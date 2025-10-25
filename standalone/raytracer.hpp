#pragma once

#include "geometry.hpp"

#include <vector>

namespace standalone {

struct RaySegment {
    std::size_t cell{0};
    double length{0.0};
};

class RayTracer {
public:
    explicit RayTracer(const RectangularMesh &mesh);

    std::vector<RaySegment> trace(const Vector3 &origin, const Vector3 &direction) const;

private:
    const RectangularMesh *mesh_{nullptr};
};

} // namespace standalone
