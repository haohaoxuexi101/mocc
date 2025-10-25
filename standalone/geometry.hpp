#pragma once

#include "material.hpp"

#include <cstddef>
#include <vector>

namespace standalone {

struct Vector3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

class RectangularMesh {
public:
    RectangularMesh() = default;
    RectangularMesh(std::vector<double> dx,
                    std::vector<double> dy,
                    std::vector<double> dz,
                    std::vector<std::size_t> material_ids,
                    MaterialLibrary materials);

    std::size_t nx() const noexcept { return nx_; }
    std::size_t ny() const noexcept { return ny_; }
    std::size_t nz() const noexcept { return nz_; }

    double x_edge(std::size_t i) const noexcept { return x_edges_[i]; }
    double y_edge(std::size_t i) const noexcept { return y_edges_[i]; }
    double z_edge(std::size_t i) const noexcept { return z_edges_[i]; }

    double width_x(std::size_t i) const noexcept { return dx_[i]; }
    double width_y(std::size_t i) const noexcept { return dy_[i]; }
    double width_z(std::size_t i) const noexcept { return dz_[i]; }

    double total_width_x() const noexcept { return x_edges_.back(); }
    double total_width_y() const noexcept { return y_edges_.back(); }
    double total_width_z() const noexcept { return z_edges_.back(); }

    std::size_t cell_count() const noexcept { return mat_ids_.size(); }

    std::size_t cell_index(std::size_t ix, std::size_t iy, std::size_t iz) const noexcept;
    const Material &material(std::size_t cell) const;
    std::size_t material_id(std::size_t cell) const noexcept { return mat_ids_[cell]; }

    Vector3 cell_center(std::size_t ix, std::size_t iy, std::size_t iz) const noexcept;

    std::size_t locate_x(double x) const;
    std::size_t locate_y(double y) const;
    std::size_t locate_z(double z) const;
    std::size_t locate_cell(const Vector3 &position) const;

    const std::vector<double> &volumes() const noexcept { return volumes_; }
    double cell_volume(std::size_t cell) const noexcept { return volumes_[cell]; }

    const MaterialLibrary &materials() const noexcept { return materials_; }

private:
    std::vector<double> dx_;
    std::vector<double> dy_;
    std::vector<double> dz_;
    std::vector<double> x_edges_;
    std::vector<double> y_edges_;
    std::vector<double> z_edges_;
    std::vector<std::size_t> mat_ids_;
    std::vector<double> volumes_;
    MaterialLibrary materials_;
    std::size_t nx_{0};
    std::size_t ny_{0};
    std::size_t nz_{0};

    void build_edges();
};

class GeometryBuilder {
public:
    GeometryBuilder(std::vector<double> dx,
                    std::vector<double> dy,
                    std::vector<double> dz,
                    std::size_t default_material);

    void fill_box(double x0,
                  double x1,
                  double y0,
                  double y1,
                  double z0,
                  double z1,
                  std::size_t material_id);

    RectangularMesh build(MaterialLibrary materials) const;

private:
    std::vector<double> dx_;
    std::vector<double> dy_;
    std::vector<double> dz_;
    std::vector<std::size_t> mat_ids_;
    std::size_t nx_{0};
    std::size_t ny_{0};
    std::size_t nz_{0};
};

std::vector<Vector3> face_centers(const RectangularMesh &mesh, int axis, bool positive_side);

} // namespace standalone
