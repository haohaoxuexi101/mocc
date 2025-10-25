#pragma once

#include "geometry.hpp"

#include <cmath>
#include <vector>

namespace standalone {

struct Direction {
    Vector3 omega{0.0, 0.0, 0.0};
    double weight{0.0};
};

inline Vector3 normalize(const Vector3 &v) {
    double mag = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (mag == 0.0) {
        return {0.0, 0.0, 0.0};
    }
    return {v.x / mag, v.y / mag, v.z / mag};
}

inline std::vector<Direction> default_directions() {
    const double four_pi = 4.0 * std::acos(-1.0);
    const double weight = four_pi / 6.0;
    return {{normalize({1.0, 0.0, 0.0}), weight},
            {normalize({-1.0, 0.0, 0.0}), weight},
            {normalize({0.0, 1.0, 0.0}), weight},
            {normalize({0.0, -1.0, 0.0}), weight},
            {normalize({0.0, 0.0, 1.0}), weight},
            {normalize({0.0, 0.0, -1.0}), weight}};
}

inline double dot(const Direction &dir, const Vector3 &normal) {
    return dir.omega.x * normal.x + dir.omega.y * normal.y + dir.omega.z * normal.z;
}

} // namespace standalone
