#include "sn_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace standalone {

namespace {
constexpr double kFourPi = 4.0 * std::acos(-1.0);
constexpr double kTiny = 1.0e-12;
}

SNSolver::SNSolver(const RectangularMesh &mesh) : mesh_(&mesh), directions_(default_directions()) {}

void SNSolver::set_directions(std::vector<Direction> directions) {
    directions_ = std::move(directions);
}

void SNSolver::set_settings(const Settings &settings) {
    settings_ = settings;
}

std::vector<double> SNSolver::solve() const {
    if (!mesh_) {
        return {};
    }

    std::vector<double> flux(mesh_->cell_count(), 0.0);
    std::vector<double> old_flux(mesh_->cell_count(), 0.0);
    std::vector<double> track_tally(mesh_->cell_count(), 0.0);

    for (std::size_t iteration = 0; iteration < settings_.max_iterations; ++iteration) {
        std::fill(track_tally.begin(), track_tally.end(), 0.0);

        for (const auto &direction : directions_) {
            double ax = std::fabs(direction.omega.x);
            double ay = std::fabs(direction.omega.y);
            double az = std::fabs(direction.omega.z);

            int axis = 0;
            double max_component = ax;
            if (ay > max_component) {
                axis = 1;
                max_component = ay;
            }
            if (az > max_component) {
                axis = 2;
                max_component = az;
            }

            if (max_component < 0.5) {
                continue;
            }

            if (axis == 0) {
                bool positive = direction.omega.x > 0.0;
                for (std::size_t iy = 0; iy < mesh_->ny(); ++iy) {
                    for (std::size_t iz = 0; iz < mesh_->nz(); ++iz) {
                        double psi = 0.0;
                        if (positive) {
                            for (std::size_t ix = 0; ix < mesh_->nx(); ++ix) {
                                std::size_t cell = mesh_->cell_index(ix, iy, iz);
                                double length = mesh_->width_x(ix);
                                const Material &mat = mesh_->material(cell);
                                double q_iso = mat.isotropic_source(old_flux[cell]) / kFourPi;
                                double sigma_t = mat.sigma_t;
                                double psi_out = psi;
                                double psi_avg = psi;
                                if (sigma_t > 0.0 && length > 0.0) {
                                    double tau = sigma_t * length;
                                    double att = std::exp(-tau);
                                    double q_div = q_iso / sigma_t;
                                    psi_out = psi * att + q_iso * (1.0 - att) / sigma_t;
                                    double term = (tau > 1.0e-8) ? (1.0 - att) / tau : (1.0 - tau * 0.5 + tau * tau / 6.0);
                                    psi_avg = q_div + (psi - q_div) * term;
                                } else {
                                    psi_out = psi + q_iso * length;
                                    psi_avg = psi + 0.5 * q_iso * length;
                                }
                                track_tally[cell] += psi_avg * direction.weight * length;
                                psi = psi_out;
                            }
                        } else {
                            for (std::size_t ix = mesh_->nx(); ix-- > 0;) {
                                std::size_t cell = mesh_->cell_index(ix, iy, iz);
                                double length = mesh_->width_x(ix);
                                const Material &mat = mesh_->material(cell);
                                double q_iso = mat.isotropic_source(old_flux[cell]) / kFourPi;
                                double sigma_t = mat.sigma_t;
                                double psi_out = psi;
                                double psi_avg = psi;
                                if (sigma_t > 0.0 && length > 0.0) {
                                    double tau = sigma_t * length;
                                    double att = std::exp(-tau);
                                    double q_div = q_iso / sigma_t;
                                    psi_out = psi * att + q_iso * (1.0 - att) / sigma_t;
                                    double term = (tau > 1.0e-8) ? (1.0 - att) / tau : (1.0 - tau * 0.5 + tau * tau / 6.0);
                                    psi_avg = q_div + (psi - q_div) * term;
                                } else {
                                    psi_out = psi + q_iso * length;
                                    psi_avg = psi + 0.5 * q_iso * length;
                                }
                                track_tally[cell] += psi_avg * direction.weight * length;
                                psi = psi_out;
                            }
                        }
                    }
                }
            } else if (axis == 1) {
                bool positive = direction.omega.y > 0.0;
                for (std::size_t ix = 0; ix < mesh_->nx(); ++ix) {
                    for (std::size_t iz = 0; iz < mesh_->nz(); ++iz) {
                        double psi = 0.0;
                        if (positive) {
                            for (std::size_t iy = 0; iy < mesh_->ny(); ++iy) {
                                std::size_t cell = mesh_->cell_index(ix, iy, iz);
                                double length = mesh_->width_y(iy);
                                const Material &mat = mesh_->material(cell);
                                double q_iso = mat.isotropic_source(old_flux[cell]) / kFourPi;
                                double sigma_t = mat.sigma_t;
                                double psi_out = psi;
                                double psi_avg = psi;
                                if (sigma_t > 0.0 && length > 0.0) {
                                    double tau = sigma_t * length;
                                    double att = std::exp(-tau);
                                    double q_div = q_iso / sigma_t;
                                    psi_out = psi * att + q_iso * (1.0 - att) / sigma_t;
                                    double term = (tau > 1.0e-8) ? (1.0 - att) / tau : (1.0 - tau * 0.5 + tau * tau / 6.0);
                                    psi_avg = q_div + (psi - q_div) * term;
                                } else {
                                    psi_out = psi + q_iso * length;
                                    psi_avg = psi + 0.5 * q_iso * length;
                                }
                                track_tally[cell] += psi_avg * direction.weight * length;
                                psi = psi_out;
                            }
                        } else {
                            for (std::size_t iy = mesh_->ny(); iy-- > 0;) {
                                std::size_t cell = mesh_->cell_index(ix, iy, iz);
                                double length = mesh_->width_y(iy);
                                const Material &mat = mesh_->material(cell);
                                double q_iso = mat.isotropic_source(old_flux[cell]) / kFourPi;
                                double sigma_t = mat.sigma_t;
                                double psi_out = psi;
                                double psi_avg = psi;
                                if (sigma_t > 0.0 && length > 0.0) {
                                    double tau = sigma_t * length;
                                    double att = std::exp(-tau);
                                    double q_div = q_iso / sigma_t;
                                    psi_out = psi * att + q_iso * (1.0 - att) / sigma_t;
                                    double term = (tau > 1.0e-8) ? (1.0 - att) / tau : (1.0 - tau * 0.5 + tau * tau / 6.0);
                                    psi_avg = q_div + (psi - q_div) * term;
                                } else {
                                    psi_out = psi + q_iso * length;
                                    psi_avg = psi + 0.5 * q_iso * length;
                                }
                                track_tally[cell] += psi_avg * direction.weight * length;
                                psi = psi_out;
                            }
                        }
                    }
                }
            } else {
                bool positive = direction.omega.z > 0.0;
                for (std::size_t ix = 0; ix < mesh_->nx(); ++ix) {
                    for (std::size_t iy = 0; iy < mesh_->ny(); ++iy) {
                        double psi = 0.0;
                        if (positive) {
                            for (std::size_t iz = 0; iz < mesh_->nz(); ++iz) {
                                std::size_t cell = mesh_->cell_index(ix, iy, iz);
                                double length = mesh_->width_z(iz);
                                const Material &mat = mesh_->material(cell);
                                double q_iso = mat.isotropic_source(old_flux[cell]) / kFourPi;
                                double sigma_t = mat.sigma_t;
                                double psi_out = psi;
                                double psi_avg = psi;
                                if (sigma_t > 0.0 && length > 0.0) {
                                    double tau = sigma_t * length;
                                    double att = std::exp(-tau);
                                    double q_div = q_iso / sigma_t;
                                    psi_out = psi * att + q_iso * (1.0 - att) / sigma_t;
                                    double term = (tau > 1.0e-8) ? (1.0 - att) / tau : (1.0 - tau * 0.5 + tau * tau / 6.0);
                                    psi_avg = q_div + (psi - q_div) * term;
                                } else {
                                    psi_out = psi + q_iso * length;
                                    psi_avg = psi + 0.5 * q_iso * length;
                                }
                                track_tally[cell] += psi_avg * direction.weight * length;
                                psi = psi_out;
                            }
                        } else {
                            for (std::size_t iz = mesh_->nz(); iz-- > 0;) {
                                std::size_t cell = mesh_->cell_index(ix, iy, iz);
                                double length = mesh_->width_z(iz);
                                const Material &mat = mesh_->material(cell);
                                double q_iso = mat.isotropic_source(old_flux[cell]) / kFourPi;
                                double sigma_t = mat.sigma_t;
                                double psi_out = psi;
                                double psi_avg = psi;
                                if (sigma_t > 0.0 && length > 0.0) {
                                    double tau = sigma_t * length;
                                    double att = std::exp(-tau);
                                    double q_div = q_iso / sigma_t;
                                    psi_out = psi * att + q_iso * (1.0 - att) / sigma_t;
                                    double term = (tau > 1.0e-8) ? (1.0 - att) / tau : (1.0 - tau * 0.5 + tau * tau / 6.0);
                                    psi_avg = q_div + (psi - q_div) * term;
                                } else {
                                    psi_out = psi + q_iso * length;
                                    psi_avg = psi + 0.5 * q_iso * length;
                                }
                                track_tally[cell] += psi_avg * direction.weight * length;
                                psi = psi_out;
                            }
                        }
                    }
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
            std::cout << "[SN] iteration " << iteration + 1 << ": max delta=" << max_delta << "\n";
        }

        old_flux = flux;

        if (max_delta < settings_.tolerance) {
            break;
        }
    }

    return flux;
}

} // namespace standalone
