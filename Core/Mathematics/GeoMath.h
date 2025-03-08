#pragma once

#include "Core/CoreExports.h"
#include "Core/Mathematics/Vector3.h"
#include <random>
#include <cmath>
#include <numbers>

namespace AeonTerra {
    namespace Math {

        struct SphericalCoord {
            float radius;
            float theta; // polar angle
            float phi;   // azimuthal angle
        };

        class AEONTERRA_API GeoMath {
        public:
            static Vector3 randomUnitVector(std::mt19937& engine) {
                std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
                while (true) {
                    Vector3 v(dist(engine), dist(engine), dist(engine));
                    if (v.lengthSquared() > 0) {
                        return v.normalized();
                    }
                }
            }

            // Implementation of the missing GreatCircleDistance function
            static float GreatCircleDistance(const SphericalCoord& a, const SphericalCoord& b) {
                // Haversine formula for great circle distance
                float sinLat1 = std::sin(a.theta);
                float sinLat2 = std::sin(b.theta);
                float cosLat1 = std::cos(a.theta);
                float cosLat2 = std::cos(b.theta);
                float cosLon = std::cos(a.phi - b.phi);

                // Calculate distance on unit sphere
                float cosDistance = sinLat1 * sinLat2 + cosLat1 * cosLat2 * cosLon;

                // Clamp to valid range to prevent numerical errors
                if (cosDistance > 1.0f) cosDistance = 1.0f;
                if (cosDistance < -1.0f) cosDistance = -1.0f;

                // Convert to actual distance using the radius
                float distance = std::acos(cosDistance) * a.radius;
                return distance;
            }
        };

    } // namespace Math
} // namespace AeonTerra