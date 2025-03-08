#pragma once
#include <cmath>
#include "Core/CoreExports.h"

namespace AeonTerra {
    namespace Math {

        class AEONTERRA_API Vector3 {
        public:
            float x, y, z;

            Vector3(float x = 0.0f, float y = 0.0f, float z = 0.0f)
                : x(x), y(y), z(z) {
            }

            float lengthSquared() const {
                return x * x + y * y + z * z;
            }

            float length() const {
                return std::sqrt(lengthSquared());
            }

            Vector3 normalized() const {
                const float len = length();
                return len > 0 ? Vector3(x / len, y / len, z / len) : *this;
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

            Vector3 operator+(const Vector3& rhs) const {
                return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
            }

            Vector3 operator-(const Vector3& rhs) const {
                return Vector3(x - rhs.x, y - rhs.y, z - rhs.z);
            }

            Vector3 cross(const Vector3& rhs) const {
                return Vector3(
                    y * rhs.z - z * rhs.y,
                    z * rhs.x - x * rhs.z,
                    x * rhs.y - y * rhs.x
                );
            }

            float dot(const Vector3& rhs) const {
                return x * rhs.x + y * rhs.y + z * rhs.z;
            }
        };

    } // namespace Math
} // namespace AeonTerra