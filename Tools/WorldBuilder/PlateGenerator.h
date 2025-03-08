#pragma once

#include "Core/CoreExports.h"
#include "Core/Mathematics/Vector3.h"
#include "Core/Mathematics/GeoMath.h"
#include <random>
#include <vector>
#include <stdexcept>

namespace AeonTerra {
namespace Tools {

struct AEONTERRA_API Plate {
    int id;
    AeonTerra::Math::Vector3 center;
    AeonTerra::Math::Vector3 velocity;
    float buoyancy;
    float thickness;
    std::vector<AeonTerra::Math::SphericalCoord> boundaries;
};

class AEONTERRA_API PlateGenerator {
public:
    explicit PlateGenerator(uint32_t seed = std::random_device{}());
    void generatePlates(int plateCount, float planetRadius);
    const std::vector<Plate>& getPlates() const;

private:
    void CreateVoronoiCells(std::vector<Plate>& plates, float sphereRadius);

    std::mt19937 m_engine;
    int m_plateCount;
    float m_planetRadius;
    std::vector<Plate> m_plates;
};

} // namespace Tools
} // namespace AeonTerra