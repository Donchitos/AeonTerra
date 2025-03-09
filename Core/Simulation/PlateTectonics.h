#pragma once
#include <vector>
#include <unordered_map>
#include "../Mathematics/Vector3.h"

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
    
    // New structure to represent boundary points (e.g., mountain ranges, rifts)
    struct BoundaryPoint {
        Math::Vector3 position;
        float height;       // Current height/elevation
        float heightChange; // Rate of change in height
    };
    
    // Collision types between plates
    enum class CollisionType {
        Transform,    // Sliding past each other
        Convergent,   // Moving toward each other
        Divergent     // Moving away from each other
    };
    
    // Structure to represent a collision between plates
    struct PlateCollision {
        int plate1Index;
        int plate2Index;
        CollisionType type;
        std::vector<Math::Vector3> points;  // Points along the collision boundary
        float relativeVelocity;
    };

    PlateTectonics();
    
    // Main simulation method - advances simulation by specified number of years
    void SimulateFrame(float deltaYears);
    
    // Load plates from the PlateGenerator tool
    void LoadPlatesFromGenerator(const Tools::PlateGenerator& generator);
    
    // Getters for simulation state
    const std::vector<Plate>& GetPlates() const { return plates; }
    float GetSimulationTime() const { return simulationTime; }
    std::vector<BoundaryPoint> GetConvergentBoundaries() const;
    
    // Set simulation parameters
    void SetMantleViscosity(float viscosity) { mantleViscosity = viscosity; }

private:
    // Calculate forces from mantle convection for each plate
    void CalculateMantleForces();
    
    // Detect all collisions between plates
    std::vector<PlateCollision> DetectAllCollisions();
    
    // Resolve identified collisions
    void ResolveCollisions(const std::vector<PlateCollision>& collisions, float deltaSeconds);
    
    // Helper method to detect collisions between plates
    bool DetectCollision(const Plate& plate1, const Plate& plate2, std::vector<Math::Vector3>& collisionPoints);
    
    // Update all convergent boundary locations
    void UpdateConvergentBoundaries(float deltaYears);
    
    // Plate data
    std::vector<Plate> plates;
    
    // Forces applied by mantle to each plate
    std::unordered_map<int, float> mantleForces;
    
    // Convergent boundary points (mountains, trenches)
    std::vector<BoundaryPoint> convergentBoundaries;
    
    // Physical parameters
    float mantleViscosity; // Pa·s
    
    // Current simulation time in years
    float simulationTime;
};

} // namespace Simulation
} // namespace AeonTerra