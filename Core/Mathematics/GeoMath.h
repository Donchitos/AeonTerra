#pragma once

namespace AeonTerra::Math {
    
struct SphericalCoord {
    float radius;
    float theta; // polar angle
    float phi;   // azimuthal angle
};

} // namespace AeonTerra::Math
#include "Core/CoreExports.h"
#include "Core/Mathematics/Vector3.h"
#include <random>
#include <numbers>

namespace AeonTerra::Math {
class AEONTERRA_API GeoMath {
public:
    static Vector3 randomUnitVector(std::mt19937& engine) {
        std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
        while(true) {
            Vector3 v(dist(engine), dist(engine), dist(engine));
            if (v.lengthSquared() > 0) {
                return v.normalized();
            }
        }
    }
};
} // namespace AeonTerra::Math