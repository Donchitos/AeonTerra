#pragma once
#include <vector>
#include <unordered_map>

// Forward declaration for PlateGenerator
namespace AeonTerra {
namespace Tools {
    class PlateGenerator;
}
}

namespace AeonTerra {
namespace Simulation {

class PlateTectonics {
public:
    struct Plate {
        int id;
        std::vector<float> boundaries;
        float angularVelocity;
        float thickness;
        float density;
    };

    PlateTectonics();
    
    // Main simulation method - advances simulation by specified number of years
    void SimulateFrame(float deltaYears);
    
    // Load plates from the PlateGenerator tool
    void LoadPlatesFromGenerator(const Tools::PlateGenerator& generator);
    
    // Getters for simulation state
    const std::vector<Plate>& GetPlates() const { return plates; }
    float GetSimulationTime() const { return simulationTime; }
    
    // Set simulation parameters
    void SetMantleViscosity(float viscosity) { mantleViscosity = viscosity; }

private:
    // Calculate forces from mantle convection for each plate
    void CalculateMantleForces();
    
    // Resolve collisions between plates
    void ResolvePlateCollisions();
    
    // Helper method to detect collisions between plates
    bool DetectCollision(const Plate& plate1, const Plate& plate2);
    
    // Plate data
    std::vector<Plate> plates;
    
    // Forces applied by mantle to each plate
    std::unordered_map<int, float> mantleForces;
    
    // Physical parameters
    float mantleViscosity; // Pa·s
    
    // Current simulation time in years
    float simulationTime;
};

} // namespace Simulation
} // namespace AeonTerra