#pragma once

#include "../../Core/CoreExports.h"
#include "../../Core/Mathematics/Vector3.h"
#include "../../Core/Mathematics/GeoMath.h"
#include "../../Core/Simulation/PlateTectonics.h"
#include <random>
#include <vector>
#include <stdexcept>

namespace AeonTerra {
namespace Tools {

// Updated Plate structure to be compatible with PlateTectonics::Plate
struct Plate {
    int id;
    Math::Vector3 center;
    Math::Vector3 velocity;
    float buoyancy;
    float thickness;
    float density;
    std::vector<Math::SphericalCoord> boundaries;
};

class PlateGenerator {
public:
    explicit PlateGenerator(uint32_t seed = std::random_device{}());
    void generatePlates(int plateCount, float planetRadius);
    const std::vector<Plate>& getPlates() const;

    // Convert to PlateTectonics::Plate format
    std::vector<Simulation::PlateTectonics::Plate> convertToTectonicPlates() const {
        std::vector<Simulation::PlateTectonics::Plate> tectonicPlates;
        for(const auto& plate : m_plates) {
            Simulation::PlateTectonics::Plate tectonicPlate;
            tectonicPlate.id = plate.id;
            
            // Convert the boundaries to float vector
            // This assumes PlateTectonics::Plate.boundaries is a vector<float>
            // and expects flattened spherical coordinates
            for(const auto& boundary : plate.boundaries) {
                tectonicPlate.boundaries.push_back(boundary.radius);
                tectonicPlate.boundaries.push_back(boundary.theta);
                tectonicPlate.boundaries.push_back(boundary.phi);
            }
            
            // Calculate angular velocity based on plate movement characteristics
            // Using m_planetRadius (class member) instead of planetRadius (which was undefined)
            tectonicPlate.angularVelocity = plate.velocity.length() / (m_planetRadius * 1000.0f); // Convert km to meters
            tectonicPlate.thickness = plate.thickness;
            tectonicPlate.density = plate.density;
            
            tectonicPlates.push_back(tectonicPlate);
        }
        return tectonicPlates;
    }

private:
    void CreateVoronoiCells(std::vector<Plate>& plates, float sphereRadius);

    std::mt19937 m_engine;
    int m_plateCount;
    float m_planetRadius;
    std::vector<Plate> m_plates;
};

} // namespace Tools
} // namespace AeonTerra