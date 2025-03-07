#pragma once
#include <vector>
#include <unordered_map>

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
    void SimulateFrame(float deltaYears);
    const std::vector<Plate>& GetPlates() const { return plates; }

private:
    void CalculateMantleForces();
    void ResolvePlateCollisions();
    
    std::vector<Plate> plates;
    std::unordered_map<int, float> mantleForces;
    float mantleViscosity = 1e21f; // Pa·s
    float simulationTime = 0.0f;
};

} // namespace Simulation
} // namespace AeonTerra