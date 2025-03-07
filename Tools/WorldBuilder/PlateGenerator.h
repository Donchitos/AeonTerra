#pragma once

#include <Core/Mathematics/Vector3.h>
#include <random>
#include <vector>
#include <stdexcept>

namespace AeonTerra {

struct Plate {
    Vector3 center;
    Vector3 velocity;
    float buoyancy;
    float thickness;
};

class PlateGenerator {
public:
    explicit PlateGenerator(uint32_t seed = std::random_device{}());
    
    void generatePlates(int plateCount, float planetRadius);
    const std::vector<Plate>& getPlates() const;

private:
    void calculatePlateBoundaries(float planetRadius);

    std::mt19937 m_engine;
    std::vector<Plate> m_plates;
};

} // namespace AeonTerra
#pragma once
#include "../Core/Simulation/PlateTectonics.h"
#include <vector>

namespace AeonTerra {
namespace Tools {

class PlateGenerator {
public:
    static std::vector<Simulation::PlateTectonics::Plate> GeneratePlates(
        int count = 12, 
        float sphereRadius = 6371.0f // Earth-like radius in km
    );
    
private:
    static void CreateVoronoiCells(std::vector<Simulation::PlateTectonics::Plate>& plates);
    static void AssignPlateProperties(std::vector<Simulation::PlateTectonics::Plate>& plates);
};

} // namespace Tools
} // namespace AeonTerra