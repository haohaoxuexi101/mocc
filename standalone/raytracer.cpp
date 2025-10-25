#include "raytracer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace standalone {

namespace {
constexpr double kEps = 1.0e-9;
constexpr double kTol = 1.0e-10;

inline double safe_ratio(double numerator, double denominator) {
    if (std::fabs(denominator) < kTol) {
        return std::numeric_limits<double>::infinity();
    }
    return numerator / denominator;
}

inline std::size_t clamp_index(std::size_t value, std::size_t max_value) {
    return (value >= max_value) ? max_value - 1 : value;
}
}

RayTracer::RayTracer(const RectangularMesh &mesh) : mesh_(&mesh) {}

std::vector<RaySegment> RayTracer::trace(const Vector3 &origin, const Vector3 &direction) const {
    if (!mesh_) {
        throw std::runtime_error("RayTracer mesh not set");
    }

    Vector3 pos = origin;
    pos.x = std::clamp(pos.x, 0.0 + kEps, mesh_->total_width_x() - kEps);
    pos.y = std::clamp(pos.y, 0.0 + kEps, mesh_->total_width_y() - kEps);
    pos.z = std::clamp(pos.z, 0.0 + kEps, mesh_->total_width_z() - kEps);

    std::size_t ix = clamp_index(mesh_->locate_x(pos.x), mesh_->nx());
    std::size_t iy = clamp_index(mesh_->locate_y(pos.y), mesh_->ny());
    std::size_t iz = clamp_index(mesh_->locate_z(pos.z), mesh_->nz());

    int step_x = (direction.x > 0.0) ? 1 : ((direction.x < 0.0) ? -1 : 0);
    int step_y = (direction.y > 0.0) ? 1 : ((direction.y < 0.0) ? -1 : 0);
    int step_z = (direction.z > 0.0) ? 1 : ((direction.z < 0.0) ? -1 : 0);

    std::vector<RaySegment> segments;
    segments.reserve(mesh_->nx() + mesh_->ny() + mesh_->nz());

    Vector3 current_pos = pos;

    while (ix < mesh_->nx() && iy < mesh_->ny() && iz < mesh_->nz()) {
        double next_x = (step_x > 0) ? mesh_->x_edge(ix + 1) : mesh_->x_edge(ix);
        double next_y = (step_y > 0) ? mesh_->y_edge(iy + 1) : mesh_->y_edge(iy);
        double next_z = (step_z > 0) ? mesh_->z_edge(iz + 1) : mesh_->z_edge(iz);

        double tx = (step_x != 0) ? safe_ratio(next_x - current_pos.x, direction.x) : std::numeric_limits<double>::infinity();
        double ty = (step_y != 0) ? safe_ratio(next_y - current_pos.y, direction.y) : std::numeric_limits<double>::infinity();
        double tz = (step_z != 0) ? safe_ratio(next_z - current_pos.z, direction.z) : std::numeric_limits<double>::infinity();

        double dt = std::min({tx, ty, tz});
        if (!std::isfinite(dt)) {
            break;
        }
        if (dt < 0.0) {
            dt = 0.0;
        }

        segments.push_back({mesh_->cell_index(ix, iy, iz), dt});

        current_pos.x += direction.x * (dt + kEps);
        current_pos.y += direction.y * (dt + kEps);
        current_pos.z += direction.z * (dt + kEps);

        bool advance_x = (tx - dt) <= kTol;
        bool advance_y = (ty - dt) <= kTol;
        bool advance_z = (tz - dt) <= kTol;

        if (advance_x) {
            if (step_x > 0) {
                ++ix;
            } else if (step_x < 0) {
                if (ix == 0) {
                    ix = mesh_->nx();
                } else {
                    --ix;
                }
            }
        }
        if (advance_y) {
            if (step_y > 0) {
                ++iy;
            } else if (step_y < 0) {
                if (iy == 0) {
                    iy = mesh_->ny();
                } else {
                    --iy;
                }
            }
        }
        if (advance_z) {
            if (step_z > 0) {
                ++iz;
            } else if (step_z < 0) {
                if (iz == 0) {
                    iz = mesh_->nz();
                } else {
                    --iz;
                }
            }
        }

        if (ix >= mesh_->nx() || iy >= mesh_->ny() || iz >= mesh_->nz()) {
            break;
        }
    }

    return segments;
}

} // namespace standalone
