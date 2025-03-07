#pragma once
#include <cmath>

namespace AeonTerra {
namespace Math {

struct SphericalCoord {
    float radius;
    float theta; // polar angle (radians)
    float phi;   // azimuthal angle (radians)
};

struct CartesianCoord {
    float x;
    float y;
    float z;
};

class GeoMath {
public:
    static CartesianCoord SphericalToCartesian(const SphericalCoord& sc) {
        return {
            sc.radius * sin(sc.theta) * cos(sc.phi),
            sc.radius * sin(sc.theta) * sin(sc.phi),
            sc.radius * cos(sc.theta)
        };
    }

    static SphericalCoord CartesianToSpherical(const CartesianCoord& cc) {
        float radius = sqrt(cc.x*cc.x + cc.y*cc.y + cc.z*cc.z);
        return {
            radius,
            acos(cc.z / radius),
            atan2(cc.y, cc.x)
        };
    }

    static float GreatCircleDistance(const SphericalCoord& a, const SphericalCoord& b) {
        return acos(
            sin(a.theta) * sin(b.theta) * cos(a.phi - b.phi) + 
            cos(a.theta) * cos(b.theta)
        );
    }
};

} // namespace Math
} // namespace AeonTerra