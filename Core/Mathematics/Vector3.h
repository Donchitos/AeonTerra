#pragma once
#include <random>
#include <cmath>
#include "Core/CoreExports.h"

namespace AeonTerra::Math {

class AEONTERRA_API Vector3 {
public:
    float x, y, z;

    Vector3(float x = 0.0f, float y = 0.0f, float z = 0.0f)
        : x(x), y(y), z(z) {}

    float lengthSquared() const {
        return x*x + y*y + z*z;
    }

    Vector3 normalized() const {
        const float len = std::sqrt(lengthSquared());
        return len > 0 ? Vector3(x/len, y/len, z/len) : *this;
    }

    Vector3 operator*(float scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    Vector3& operator+=(const Vector3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    static Vector3 randomUnitVector(std::mt19937& engine);
};

} // namespace AeonTerra::Math