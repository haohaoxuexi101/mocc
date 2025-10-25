#include "geometry.hpp"
#include "material.hpp"
#include "moc_solver.hpp"
#include "sn_solver.hpp"
#include "transport_common.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

using namespace standalone;

int main() {
    MaterialLibrary library;
    std::size_t moderator = library.add_material("moderator", Material(0.3, 0.2, 0.0, 1.0, 0.1));
    std::size_t fuel = library.add_material("fuel", Material(0.8, 0.6, 0.0, 1.0, 0.5));

    std::vector<double> spacing(4, 1.0);
    std::vector<double> axial(3, 1.5);

    GeometryBuilder builder(spacing, spacing, axial, moderator);
    builder.fill_box(1.0, 3.0, 1.0, 3.0, 0.0, axial.back(), fuel);

    RectangularMesh mesh = builder.build(library);

    MOCSolver moc(mesh);
    MOCSolver::Settings moc_settings;
    moc_settings.max_iterations = 400;
    moc_settings.tolerance = 1.0e-5;
    moc.set_settings(moc_settings);
    auto moc_flux = moc.solve();

    SNSolver sn(mesh);
    SNSolver::Settings sn_settings;
    sn_settings.max_iterations = 400;
    sn_settings.tolerance = 1.0e-5;
    sn.set_settings(sn_settings);
    auto sn_flux = sn.solve();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Cell index | MoC flux | Sn flux" << std::endl;
    for (std::size_t idx = 0; idx < mesh.cell_count(); ++idx) {
        std::cout << std::setw(10) << idx << " | " << std::setw(9) << moc_flux[idx] << " | "
                  << std::setw(8) << sn_flux[idx] << '\n';
    }

    return 0;
}
