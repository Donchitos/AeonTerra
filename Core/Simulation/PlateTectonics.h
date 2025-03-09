#pragma once
#include <vector>
#include <unordered_map>
#include <memory>
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
        bool isOceanic;
        float age;
        std::vector<int> neighborIds;
    };
    
    // Enhanced boundary point with more properties
    struct BoundaryPoint {
        Math::Vector3 position;
        float height;       // Current height/elevation
        float heightChange; // Rate of change in height
        int plate1Id;       // First plate ID
        int plate2Id;       // Second plate ID (-1 if not a boundary)
        float divergenceRate; // Rate of movement (+ for divergent, - for convergent)
        float slipRate;     // Rate of transform motion
        float age;          // Age of the boundary
    };
    
    // Mantle cell structure for improved convection modeling
    struct MantelCell {
        Math::Vector3 center;
        Math::Vector3 flowDirection;
        float strength;
        float temperature;
        float depth;
    };
    
    // Force vector with direction
    struct ForceVector {
        float magnitude;
        Math::Vector3 direction;
    };
    
    // Enhanced collision types
    enum class CollisionType {
        Transform,    // Sliding past each other
        Convergent,   // Moving toward each other
        Divergent,    // Moving away from each other
        ObliqueConvergent, // Convergent with lateral component
        ObliqueDivergent   // Divergent with lateral component
    };
    
    // Structure to represent a collision between plates
    struct PlateCollision {
        int plate1Index;
        int plate2Index;
        CollisionType type;
        std::vector<Math::Vector3> points;  // Points along the collision boundary
        float relativeVelocity;
        float convergenceAngle;  // Angle of convergence/divergence
        bool oceanicSubduction;  // Whether oceanic plate is subducting
        float slipRate;          // Rate of transform motion
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
    const std::vector<MantelCell>& GetMantleCells() const { return mantleCells; }
    
    // Set simulation parameters
    void SetMantleViscosity(float viscosity) { mantleViscosity = viscosity; }
    void SetSeafloorSpreadingRate(float rate) { seafloorSpreadingRate = rate; }
    void SetSubductionRate(float rate) { subductionRate = rate; }
    void SetContinentalCollisionRate(float rate) { continentalCollisionRate = rate; }

private:
    // New improved methods for realistic simulation
    void CalculateImprovedMantleForces();
    void UpdatePlateAngularMomentum(float deltaSeconds);
    void ClassifyPlateBoundaries();
    
    // Enhanced boundary processing
    void ProcessSubductionZones(float deltaSeconds);
    void ProcessSpreadingRidges(float deltaSeconds);
    void ProcessTransformFaults(float deltaSeconds);
    void ProcessContinentalCollisions(float deltaSeconds);
    
    // Helper methods
    std::vector<MantelCell> CreateEarthLikeConvectionPattern();
    std::vector<Math::Vector3> SamplePlatePoints(const Plate& plate, int sampleCount);
    float CalculateCellInfluence(float distance, const MantelCell& cell);
    
    // Original methods (kept for compatibility)
    void CalculateMantleForces();
    std::vector<PlateCollision> DetectAllCollisions();
    void ResolveCollisions(const std::vector<PlateCollision>& collisions, float deltaSeconds);
    bool DetectCollision(const Plate& plate1, const Plate& plate2, std::vector<Math::Vector3>& collisionPoints);
    void UpdateConvergentBoundaries(float deltaYears);
    
    // Plate data
    std::vector<Plate> plates;
    
    // Forces applied by mantle to each plate
    std::unordered_map<int, float> mantleForces;
    std::unordered_map<int, ForceVector> mantleForcesWithDirection;
    
    // Convergent boundary points (mountains, trenches)
    std::vector<BoundaryPoint> convergentBoundaries;
    
    // New: Classified boundary points
    std::vector<BoundaryPoint> divergentBoundaries;
    std::vector<BoundaryPoint> transformBoundaries;
    
    // Mantle convection cells
    std::vector<MantelCell> mantleCells;
    
    // Physical parameters
    float mantleViscosity; // Pa·s
    float seafloorSpreadingRate; // mm/year
    float subductionRate; // mm/year
    float continentalCollisionRate; // mm/year
    
    // Current simulation time in years
    float simulationTime;
};

} // namespace Simulation
} // namespace AeonTerra