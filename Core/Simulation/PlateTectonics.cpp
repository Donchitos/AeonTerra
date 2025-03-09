#include "PlateTectonics.h"
#include "../../Tools/WorldBuilder/PlateGenerator.h"
#include <cmath>
#include <algorithm>
#include <random>

namespace AeonTerra {
namespace Simulation {

PlateTectonics::PlateTectonics() 
    : simulationTime(0.0f), mantleViscosity(1e21f) {
    // Initialize with empty plates array
}

void PlateTectonics::SimulateFrame(float deltaYears) {
    // Convert years to seconds for physics calculations
    const float deltaSeconds = deltaYears * 365.25f * 24.0f * 60.0f * 60.0f;
    
    // Calculate forces from mantle convection
    CalculateMantleForces();
    
    // Apply forces to plates and update velocities
    for (auto& plate : plates) {
        // Apply mantle force to angular velocity
        float mantleForce = mantleForces[plate.id];
        plate.angularVelocity += mantleForce * deltaSeconds / (plate.thickness * plate.density);
        
        // Apply viscous drag from mantle (simplified model)
        float dragCoefficient = 0.05f;
        plate.angularVelocity *= (1.0f - dragCoefficient * deltaSeconds);
    }
    
    // Resolve collisions between plates
    ResolvePlateCollisions();
    
    // Update simulation time
    simulationTime += deltaYears;
}

void PlateTectonics::CalculateMantleForces() {
    // Clear previous forces
    mantleForces.clear();
    
    // Random number generation for stochastic mantle convection
    std::mt19937 rng(static_cast<unsigned int>(simulationTime * 1000.0f));
    std::normal_distribution<float> forceDist(0.0f, 1e15f); // Mean force is zero, std dev based on typical mantle forces
    
    // Calculate new mantle forces for each plate
    for (const auto& plate : plates) {
        // Generate a force based on current simulation time and plate properties
        float baseForce = forceDist(rng);
        
        // Adjust force based on plate thickness (thicker plates receive less force)
        float thicknessModifier = 100.0f / (plate.thickness + 1.0f);
        
        // Density affects buoyancy in mantle
        float densityModifier = 3000.0f / plate.density;
        
        // Calculate final force
        float totalForce = baseForce * thicknessModifier * densityModifier;
        
        // Store force
        mantleForces[plate.id] = totalForce;
    }
}

void PlateTectonics::ResolvePlateCollisions() {
    // Detect plate boundary collisions using relative velocities
    for (auto& plate1 : plates) {
        for (auto& plate2 : plates) {
            if (plate1.id >= plate2.id) continue; // Skip self-collision and duplicate pairs
            
            // Check if plates are potentially colliding (simplified check)
            bool collision = DetectCollision(plate1, plate2);
            
            if (collision) {
                // Determine collision type based on relative motion
                float relativeVelocity = plate1.angularVelocity - plate2.angularVelocity;
                
                if (std::abs(relativeVelocity) < 0.00001f) {
                    // Transform boundary (plates sliding past each other)
                    // Create fault lines but minimal vertical displacement
                    // In this simplified model, we just slow both plates slightly
                    plate1.angularVelocity *= 0.95f;
                    plate2.angularVelocity *= 0.95f;
                }
                else if (relativeVelocity > 0) {
                    // Convergent boundary (plates moving toward each other)
                    if (plate1.density > plate2.density) {
                        // Subduction - denser plate (plate1) slides under less dense (plate2)
                        plate1.thickness += plate2.thickness * 0.05f;
                        plate2.angularVelocity *= 0.9f;
                        // In a more complex model, we would also adjust the boundaries
                    } else {
                        // Mountain building - uplift both plates
                        plate1.thickness *= 1.05f;
                        plate2.thickness *= 1.05f;
                        // Collision slows both plates
                        plate1.angularVelocity *= 0.8f;
                        plate2.angularVelocity *= 0.8f;
                    }
                }
                else {
                    // Divergent boundary (plates moving away from each other)
                    // New crust forms at the boundary
                    // Simplified: thin both plates slightly as material rises from mantle
                    plate1.thickness *= 0.99f;
                    plate2.thickness *= 0.99f;
                    // Slightly increase velocity due to reduced friction
                    plate1.angularVelocity *= 1.01f;
                    plate2.angularVelocity *= 1.01f;
                }
            }
        }
    }
}

bool PlateTectonics::DetectCollision(const Plate& plate1, const Plate& plate2) {
    // This is a simplified collision detection algorithm
    // In a real simulation, we would use the boundaries and check for intersections
    
    // For now, we'll use a simple probability model based on angular velocities
    float relativeVelocity = std::abs(plate1.angularVelocity - plate2.angularVelocity);
    float collisionProbability = relativeVelocity * 10.0f; // Higher relative velocity = higher collision chance
    
    // Clamp probability
    collisionProbability = std::min(collisionProbability, 0.5f);
    
    // Random check for collision
    return (static_cast<float>(rand()) / RAND_MAX) < collisionProbability;
}

void PlateTectonics::LoadPlatesFromGenerator(const Tools::PlateGenerator& generator) {
    // Convert generator plates to tectonic plates
    plates = generator.convertToTectonicPlates();
    
    // Initialize mantle forces for all plates
    for (const auto& plate : plates) {
        mantleForces[plate.id] = 0.0f;
    }
}

} // namespace Simulation
} // namespace AeonTerra