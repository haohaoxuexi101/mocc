#pragma once

#include "geometry.hpp"
#include "transport_common.hpp"

#include <cstddef>
#include <vector>

namespace standalone {

class SNSolver {
public:
    struct Settings {
        std::size_t max_iterations{200};
        double tolerance{1.0e-6};
        bool verbose{false};
    };

    explicit SNSolver(const RectangularMesh &mesh);

    void set_directions(std::vector<Direction> directions);
    void set_settings(const Settings &settings);

    std::vector<double> solve() const;

private:
    const RectangularMesh *mesh_{nullptr};
    std::vector<Direction> directions_;
    Settings settings_;
};

} // namespace standalone
