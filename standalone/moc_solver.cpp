#include "moc_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace standalone {

namespace {
constexpr double kFourPi = 4.0 * std::acos(-1.0);
constexpr double kTiny = 1.0e-12;
}

MOCSolver::MOCSolver(const RectangularMesh &mesh) : mesh_(&mesh), directions_(default_directions()) {}

void MOCSolver::set_directions(std::vector<Direction> directions) {
    directions_ = std::move(directions);
}

void MOCSolver::set_settings(const Settings &settings) {
    settings_ = settings;
}

std::vector<Vector3> MOCSolver::entry_points(const Direction &direction) const {
    double ax = std::fabs(direction.omega.x);
    double ay = std::fabs(direction.omega.y);
    double az = std::fabs(direction.omega.z);

    int axis = 0;
    double best = ax;
    if (ay > best) {
        axis = 1;
        best = ay;
    }
    if (az > best) {
        axis = 2;
    }

    bool positive_side = false;
    if (axis == 0) {
        positive_side = direction.omega.x < 0.0;
    } else if (axis == 1) {
        positive_side = direction.omega.y < 0.0;
    } else {
        positive_side = direction.omega.z < 0.0;
    }

    return face_centers(*mesh_, axis, positive_side);
}

std::vector<double> MOCSolver::solve() const {
    if (!mesh_) {
        return {};
    }

    std::vector<double> flux(mesh_->cell_count(), 0.0);
    std::vector<double> old_flux(mesh_->cell_count(), 0.0);
    std::vector<double> track_tally(mesh_->cell_count(), 0.0);

    RayTracer tracer(*mesh_);

    for (std::size_t iteration = 0; iteration < settings_.max_iterations; ++iteration) {
        std::fill(track_tally.begin(), track_tally.end(), 0.0);

        for (const auto &direction : directions_) {
            auto starts = entry_points(direction);
            for (const auto &start : starts) {
                double psi_in = 0.0;
                auto segments = tracer.trace(start, direction.omega);
                for (const auto &segment : segments) {
                    const Material &mat = mesh_->material(segment.cell);
                    double sigma_t = mat.sigma_t;
                    double q_iso = mat.isotropic_source(old_flux[segment.cell]) / kFourPi;
                    double psi_out = psi_in;
                    double psi_avg = psi_in;

                    if (segment.length > 0.0) {
                        if (sigma_t > 0.0) {
                            double tau = sigma_t * segment.length;
                            double att = std::exp(-tau);
                            double q_div = q_iso / sigma_t;
                            psi_out = psi_in * att + q_iso * (1.0 - att) / sigma_t;
                            double term = (tau > 1.0e-8) ? (1.0 - att) / tau : (1.0 - tau * 0.5 + tau * tau / 6.0);
                            psi_avg = q_div + (psi_in - q_div) * term;
                        } else {
                            psi_out = psi_in + q_iso * segment.length;
                            psi_avg = psi_in + 0.5 * q_iso * segment.length;
                        }
                    }

                    track_tally[segment.cell] += psi_avg * direction.weight * segment.length;
                    psi_in = psi_out;
                }
            }
        }

        double max_delta = 0.0;
        for (std::size_t cell = 0; cell < flux.size(); ++cell) {
            double volume = mesh_->cell_volume(cell);
            double new_flux = (volume > kTiny) ? track_tally[cell] / volume : 0.0;
            double delta = std::fabs(new_flux - old_flux[cell]);
            if (old_flux[cell] > kTiny) {
                delta /= old_flux[cell];
            }
            max_delta = std::max(max_delta, delta);
            flux[cell] = new_flux;
        }

        if (settings_.verbose) {
            std::cout << "[MOC] iteration " << iteration + 1 << ": max delta=" << max_delta << "\n";
        }

        old_flux = flux;

        if (max_delta < settings_.tolerance) {
            break;
        }
    }

    return flux;
}

} // namespace standalone
