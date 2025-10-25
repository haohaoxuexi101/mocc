#include "geometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace standalone {

namespace {
std::vector<double> build_edges_from_widths(const std::vector<double> &widths) {
    std::vector<double> edges(widths.size() + 1, 0.0);
    for (std::size_t i = 0; i < widths.size(); ++i) {
        edges[i + 1] = edges[i] + widths[i];
    }
    return edges;
}
}

RectangularMesh::RectangularMesh(std::vector<double> dx,
                                 std::vector<double> dy,
                                 std::vector<double> dz,
                                 std::vector<std::size_t> material_ids,
                                 MaterialLibrary materials)
    : dx_(std::move(dx)),
      dy_(std::move(dy)),
      dz_(std::move(dz)),
      mat_ids_(std::move(material_ids)),
      materials_(std::move(materials)) {
    nx_ = dx_.size();
    ny_ = dy_.size();
    nz_ = dz_.size();
    if (nx_ == 0 || ny_ == 0 || nz_ == 0) {
        throw std::runtime_error("RectangularMesh requires positive dimensions");
    }
    if (mat_ids_.size() != nx_ * ny_ * nz_) {
        throw std::runtime_error("material id array has incorrect size");
    }
    build_edges();
    volumes_.resize(mat_ids_.size());
    for (std::size_t ix = 0; ix < nx_; ++ix) {
        for (std::size_t iy = 0; iy < ny_; ++iy) {
            for (std::size_t iz = 0; iz < nz_; ++iz) {
                std::size_t idx = cell_index(ix, iy, iz);
                volumes_[idx] = dx_[ix] * dy_[iy] * dz_[iz];
            }
        }
    }
}

void RectangularMesh::build_edges() {
    x_edges_ = build_edges_from_widths(dx_);
    y_edges_ = build_edges_from_widths(dy_);
    z_edges_ = build_edges_from_widths(dz_);
}

std::size_t RectangularMesh::cell_index(std::size_t ix, std::size_t iy, std::size_t iz) const noexcept {
    return (ix * ny_ + iy) * nz_ + iz;
}

const Material &RectangularMesh::material(std::size_t cell) const {
    return materials_.material(mat_ids_[cell]);
}

Vector3 RectangularMesh::cell_center(std::size_t ix, std::size_t iy, std::size_t iz) const noexcept {
    Vector3 p;
    p.x = 0.5 * (x_edges_[ix] + x_edges_[ix + 1]);
    p.y = 0.5 * (y_edges_[iy] + y_edges_[iy + 1]);
    p.z = 0.5 * (z_edges_[iz] + z_edges_[iz + 1]);
    return p;
}

std::size_t RectangularMesh::locate_x(double x) const {
    if (x < x_edges_.front() || x > x_edges_.back()) {
        throw std::out_of_range("x coordinate outside mesh");
    }
    if (x == x_edges_.back()) {
        return nx_ - 1;
    }
    auto iter = std::upper_bound(x_edges_.begin(), x_edges_.end(), x);
    return static_cast<std::size_t>((iter - x_edges_.begin()) - 1);
}

std::size_t RectangularMesh::locate_y(double y) const {
    if (y < y_edges_.front() || y > y_edges_.back()) {
        throw std::out_of_range("y coordinate outside mesh");
    }
    if (y == y_edges_.back()) {
        return ny_ - 1;
    }
    auto iter = std::upper_bound(y_edges_.begin(), y_edges_.end(), y);
    return static_cast<std::size_t>((iter - y_edges_.begin()) - 1);
}

std::size_t RectangularMesh::locate_z(double z) const {
    if (z < z_edges_.front() || z > z_edges_.back()) {
        throw std::out_of_range("z coordinate outside mesh");
    }
    if (z == z_edges_.back()) {
        return nz_ - 1;
    }
    auto iter = std::upper_bound(z_edges_.begin(), z_edges_.end(), z);
    return static_cast<std::size_t>((iter - z_edges_.begin()) - 1);
}

std::size_t RectangularMesh::locate_cell(const Vector3 &position) const {
    std::size_t ix = locate_x(position.x);
    std::size_t iy = locate_y(position.y);
    std::size_t iz = locate_z(position.z);
    return cell_index(ix, iy, iz);
}

GeometryBuilder::GeometryBuilder(std::vector<double> dx,
                                 std::vector<double> dy,
                                 std::vector<double> dz,
                                 std::size_t default_material)
    : dx_(std::move(dx)), dy_(std::move(dy)), dz_(std::move(dz)) {
    nx_ = dx_.size();
    ny_ = dy_.size();
    nz_ = dz_.size();
    if (nx_ == 0 || ny_ == 0 || nz_ == 0) {
        throw std::runtime_error("GeometryBuilder requires non-empty widths");
    }
    mat_ids_.assign(nx_ * ny_ * nz_, default_material);
}

void GeometryBuilder::fill_box(double x0,
                               double x1,
                               double y0,
                               double y1,
                               double z0,
                               double z1,
                               std::size_t material_id) {
    std::vector<double> x_edges = build_edges_from_widths(dx_);
    std::vector<double> y_edges = build_edges_from_widths(dy_);
    std::vector<double> z_edges = build_edges_from_widths(dz_);

    for (std::size_t ix = 0; ix < nx_; ++ix) {
        double cx = 0.5 * (x_edges[ix] + x_edges[ix + 1]);
        if (cx < x0 || cx > x1) {
            continue;
        }
        for (std::size_t iy = 0; iy < ny_; ++iy) {
            double cy = 0.5 * (y_edges[iy] + y_edges[iy + 1]);
            if (cy < y0 || cy > y1) {
                continue;
            }
            for (std::size_t iz = 0; iz < nz_; ++iz) {
                double cz = 0.5 * (z_edges[iz] + z_edges[iz + 1]);
                if (cz < z0 || cz > z1) {
                    continue;
                }
                std::size_t idx = (ix * ny_ + iy) * nz_ + iz;
                mat_ids_[idx] = material_id;
            }
        }
    }
}

RectangularMesh GeometryBuilder::build(MaterialLibrary materials) const {
    return RectangularMesh(dx_, dy_, dz_, mat_ids_, std::move(materials));
}

std::vector<Vector3> face_centers(const RectangularMesh &mesh, int axis, bool positive_side) {
    std::vector<Vector3> centers;
    if (axis == 0) {
        centers.reserve(mesh.ny() * mesh.nz());
        double x = positive_side ? mesh.total_width_x() : 0.0;
        if (!positive_side) {
            x += 1e-9;
        } else {
            x -= 1e-9;
        }
        for (std::size_t iy = 0; iy < mesh.ny(); ++iy) {
            double y = 0.5 * (mesh.y_edge(iy) + mesh.y_edge(iy + 1));
            for (std::size_t iz = 0; iz < mesh.nz(); ++iz) {
                double z = 0.5 * (mesh.z_edge(iz) + mesh.z_edge(iz + 1));
                centers.push_back({x, y, z});
            }
        }
    } else if (axis == 1) {
        centers.reserve(mesh.nx() * mesh.nz());
        double y = positive_side ? mesh.total_width_y() : 0.0;
        if (!positive_side) {
            y += 1e-9;
        } else {
            y -= 1e-9;
        }
        for (std::size_t ix = 0; ix < mesh.nx(); ++ix) {
            double x = 0.5 * (mesh.x_edge(ix) + mesh.x_edge(ix + 1));
            for (std::size_t iz = 0; iz < mesh.nz(); ++iz) {
                double z = 0.5 * (mesh.z_edge(iz) + mesh.z_edge(iz + 1));
                centers.push_back({x, y, z});
            }
        }
    } else {
        centers.reserve(mesh.nx() * mesh.ny());
        double z = positive_side ? mesh.total_width_z() : 0.0;
        if (!positive_side) {
            z += 1e-9;
        } else {
            z -= 1e-9;
        }
        for (std::size_t ix = 0; ix < mesh.nx(); ++ix) {
            double x = 0.5 * (mesh.x_edge(ix) + mesh.x_edge(ix + 1));
            for (std::size_t iy = 0; iy < mesh.ny(); ++iy) {
                double y = 0.5 * (mesh.y_edge(iy) + mesh.y_edge(iy + 1));
                centers.push_back({x, y, z});
            }
        }
    }
    return centers;
}

} // namespace standalone
