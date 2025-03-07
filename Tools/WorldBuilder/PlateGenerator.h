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