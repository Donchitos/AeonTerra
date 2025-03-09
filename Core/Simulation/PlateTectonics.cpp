#include "PlateTectonics.h"
#include "../../Tools/WorldBuilder/PlateGenerator.h"
#include "../Mathematics/GeoMath.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <unordered_set>

// Include math.h for constants
#include <math.h>

// Fallback definition if still not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace AeonTerra {
namespace Simulation {

PlateTectonics::PlateTectonics() 
    : simulationTime(0.0f), mantleViscosity(1e21f) {
    // Initialize with empty plates array
}

// MODIFY the existing SimulateFrame method - don't add a new one
void PlateTectonics::SimulateFrame(float deltaYears) {
    // Convert years to seconds for physics calculations
    const float deltaSeconds = deltaYears * 365.25f * 24.0f * 60.0f * 60.0f;
    
    // Classify boundaries into convergent, divergent, and transform
    ClassifyPlateBoundaries();
    
    // Calculate forces from mantle convection using improved model
    CalculateImprovedMantleForces();
    
    // Process different boundary types
    ProcessSubductionZones(deltaSeconds);
    ProcessSpreadingRidges(deltaSeconds);
    ProcessTransformFaults(deltaSeconds);
    ProcessContinentalCollisions(deltaSeconds);
    
    // Update plate angular momentum
    UpdatePlateAngularMomentum(deltaSeconds);
    
    // Age plates
    for (auto& plate : plates) {
        plate.age += deltaYears;
    }
    
    // Update simulation time
    simulationTime += deltaYears;
}

void PlateTectonics::CalculateMantleForces() {
    // Clear previous forces
    mantleForces.clear();
    
    // Improved mantle convection model using simulated heat flow
    // We'll create a basic convection cells pattern across the sphere
    const int numConvectionCells = 12; // Approx. matching Earth's major convection cells
    
    for (const auto& plate : plates) {
        float totalForce = 0.0f;
        
        // For each plate, we'll check how it's positioned relative to convection cells
        for (size_t i = 0; i < plate.boundaries.size(); i += 3) {
            if (i + 2 >= plate.boundaries.size()) break;
            
            // Get the spherical coordinate at this boundary point
            Math::SphericalCoord point{
                plate.boundaries[i],
                plate.boundaries[i+1],
                plate.boundaries[i+2]
            };
            
            // Convert to Cartesian for easier math
            Math::Vector3 cartesian(
                point.radius * sin(point.theta) * cos(point.phi),
                point.radius * sin(point.theta) * sin(point.phi),
                point.radius * cos(point.theta)
            );
            
            // Calculate influence from each convection cell
            for (int cell = 0; cell < numConvectionCells; cell++) {
                // Simplistic model: convection cells are evenly spaced on sphere
                float cellTheta = static_cast<float>(M_PI) * (0.25f + (float)cell / numConvectionCells);
                float cellPhi = 2.0f * static_cast<float>(M_PI) * (float)cell / numConvectionCells;
                
                Math::Vector3 cellCenter(
                    sin(cellTheta) * cos(cellPhi),
                    sin(cellTheta) * sin(cellPhi),
                    cos(cellTheta)
                );
                
                // Distance to convection cell center (normalized)
                float distance = (cartesian - cellCenter).length();
                
                // Convection force falls off with distance, but reverses at cell edges
                float forceMagnitude = sin(distance * static_cast<float>(M_PI)) * 1e15f; // Scale to appropriate magnitude
                
                // Add to total force on this plate
                totalForce += forceMagnitude;
            }
        }
        
        // Average the force by number of boundary points sampled
        if (!plate.boundaries.empty()) {
            totalForce /= (plate.boundaries.size() / 3);
        }
        
        // Apply plate property modifiers
        // Thicker plates receive less force
        float thicknessModifier = 100.0f / (plate.thickness + 1.0f);
        
        // Density affects buoyancy in mantle
        float densityModifier = 3000.0f / plate.density;
        
        // Calculate final force
        float appliedForce = totalForce * thicknessModifier * densityModifier;
        
        // Store force
        mantleForces[plate.id] = appliedForce;
    }
}

std::vector<PlateTectonics::PlateCollision> PlateTectonics::DetectAllCollisions() {
    std::vector<PlateCollision> collisions;
    
    // Check each pair of plates for collisions
    for (size_t i = 0; i < plates.size(); i++) {
        for (size_t j = i + 1; j < plates.size(); j++) {
            const Plate& plate1 = plates[i];
            const Plate& plate2 = plates[j];
            
            // Detect collision between these plates
            std::vector<Math::Vector3> collisionPoints;
            float relativeVelocity = plate1.angularVelocity - plate2.angularVelocity;
            
            if (DetectCollision(plate1, plate2, collisionPoints) && !collisionPoints.empty()) {
                // Determine collision type based on relative motion
                CollisionType collisionType;
                
                if (std::abs(relativeVelocity) < 0.00001f) {
                    collisionType = CollisionType::Transform;
                }
                else if (relativeVelocity > 0) {
                    collisionType = CollisionType::Convergent;
                }
                else {
                    collisionType = CollisionType::Divergent;
                }
                
                // Add to collisions list
                PlateCollision collision;
                collision.plate1Index = static_cast<int>(i);
                collision.plate2Index = static_cast<int>(j);
                collision.type = collisionType;
                collision.points = collisionPoints;
                collision.relativeVelocity = relativeVelocity;
                
                collisions.push_back(collision);
            }
        }
    }
    
    return collisions;
}

void PlateTectonics::ResolveCollisions(const std::vector<PlateCollision>& collisions, float deltaSeconds) {
    // Remember which plates have been modified to avoid double-modification
    std::unordered_set<int> modifiedPlates;
    
    // Store boundary modifications for later processing
    struct BoundaryModification {
        int plateIndex;
        Math::Vector3 point;
        float heightChange;
    };
    std::vector<BoundaryModification> modifications;
    
    // Process each collision
    for (const auto& collision : collisions) {
        // Ensure indices are valid
        if (collision.plate1Index < 0 || collision.plate1Index >= static_cast<int>(plates.size()) ||
            collision.plate2Index < 0 || collision.plate2Index >= static_cast<int>(plates.size())) {
            continue; // Skip this collision if indices are out of bounds
        }
        
        Plate& plate1 = plates[static_cast<size_t>(collision.plate1Index)];
        Plate& plate2 = plates[static_cast<size_t>(collision.plate2Index)];
        
        switch (collision.type) {
            case CollisionType::Transform: {
                // Transform boundary (plates sliding past each other)
                // Create fault lines but minimal vertical displacement
                if (modifiedPlates.find(plate1.id) == modifiedPlates.end()) {
                    plate1.angularVelocity *= 0.95f;
                    modifiedPlates.insert(plate1.id);
                }
                
                if (modifiedPlates.find(plate2.id) == modifiedPlates.end()) {
                    plate2.angularVelocity *= 0.95f;
                    modifiedPlates.insert(plate2.id);
                }
                
                // Record small changes along transform boundary
                for (const auto& point : collision.points) {
                    modifications.push_back({static_cast<int>(collision.plate1Index), point, 100.0f * deltaSeconds});
                    modifications.push_back({static_cast<int>(collision.plate2Index), point, 100.0f * deltaSeconds});
                }
                break;
            }
            
            case CollisionType::Convergent: {
                // Convergent boundary (plates moving toward each other)
                if (plate1.density > plate2.density) {
                    // Subduction - denser plate (plate1) slides under less dense (plate2)
                    // Increase thickness of subducting plate
                    if (modifiedPlates.find(plate1.id) == modifiedPlates.end()) {
                        plate1.thickness += plate2.thickness * 0.05f * deltaSeconds;
                        modifiedPlates.insert(plate1.id);
                    }
                    
                    // Slow down overriding plate
                    if (modifiedPlates.find(plate2.id) == modifiedPlates.end()) {
                        plate2.angularVelocity *= 0.9f;
                        modifiedPlates.insert(plate2.id);
                    }
                    
                    // Record trench and volcanic arc formation
                    for (const auto& point : collision.points) {
                        // Oceanic trench at the boundary
                        modifications.push_back({static_cast<int>(collision.plate1Index), point, -2000.0f * deltaSeconds});
                        // Volcanic arc on the overriding plate
                        Math::Vector3 arcPoint = point + (point * 0.05f); // Slightly inward
                        modifications.push_back({static_cast<int>(collision.plate2Index), arcPoint, 1000.0f * deltaSeconds});
                    }
                } 
                else {
                    // Continental collision - mountain building
                    if (modifiedPlates.find(plate1.id) == modifiedPlates.end()) {
                        plate1.thickness *= (1.0f + 0.05f * deltaSeconds);
                        plate1.angularVelocity *= 0.8f;
                        modifiedPlates.insert(plate1.id);
                    }
                    
                    if (modifiedPlates.find(plate2.id) == modifiedPlates.end()) {
                        plate2.thickness *= (1.0f + 0.05f * deltaSeconds);
                        plate2.angularVelocity *= 0.8f;
                        modifiedPlates.insert(plate2.id);
                    }
                    
                    // Record mountain formation
                    for (const auto& point : collision.points) {
                        float heightIncrease = 1000.0f * deltaSeconds * std::abs(collision.relativeVelocity) * 10000.0f;
                        modifications.push_back({static_cast<int>(collision.plate1Index), point, heightIncrease});
                        modifications.push_back({static_cast<int>(collision.plate2Index), point, heightIncrease});
                    }
                }
                break;
            }
            
            case CollisionType::Divergent: {
                // Divergent boundary (plates moving away from each other)
                // New crust forms at the boundary
                if (modifiedPlates.find(plate1.id) == modifiedPlates.end()) {
                    plate1.thickness *= (1.0f - 0.01f * deltaSeconds);
                    plate1.angularVelocity *= (1.0f + 0.01f * deltaSeconds);
                    modifiedPlates.insert(plate1.id);
                }
                
                if (modifiedPlates.find(plate2.id) == modifiedPlates.end()) {
                    plate2.thickness *= (1.0f - 0.01f * deltaSeconds);
                    plate2.angularVelocity *= (1.0f + 0.01f * deltaSeconds);
                    modifiedPlates.insert(plate2.id);
                }
                
                // Record rift formation
                for (const auto& point : collision.points) {
                    modifications.push_back({static_cast<int>(collision.plate1Index), point, -500.0f * deltaSeconds});
                    modifications.push_back({static_cast<int>(collision.plate2Index), point, -500.0f * deltaSeconds});
                }
                break;
            }
        }
    }
    
    // Apply all boundary modifications
    // In a full implementation, this would update a height field or other data structure
    // For now, we'll store the modifications in the convergentBoundaries member
    for (const auto& mod : modifications) {
        convergentBoundaries.push_back({
            mod.point,
            mod.heightChange
        });
    }
}

bool PlateTectonics::DetectCollision(const Plate& plate1, const Plate& plate2, 
                                     std::vector<Math::Vector3>& collisionPoints) {
    // This is a more sophisticated collision detection algorithm
    // that identifies specific collision points between plate boundaries
    collisionPoints.clear();
    
    // Early exit if either plate has no boundaries
    if (plate1.boundaries.size() < 3 || plate2.boundaries.size() < 3)
        return false;
    
    // Assume these are in groups of 3 (radius, theta, phi)
    const float collisionThreshold = 0.05f; // 5% of planet radius as threshold
    
    // Check each point on plate1 against plate2 boundary
    for (size_t i = 0; i < plate1.boundaries.size(); i += 3) {
        if (i + 2 >= plate1.boundaries.size()) break;
        
        Math::SphericalCoord p1Point{
            plate1.boundaries[i],
            plate1.boundaries[i+1],
            plate1.boundaries[i+2]
        };
        
        // Convert to cartesian for easier distance calculation
        Math::Vector3 p1CartPoint(
            p1Point.radius * sin(p1Point.theta) * cos(p1Point.phi),
            p1Point.radius * sin(p1Point.theta) * sin(p1Point.phi),
            p1Point.radius * cos(p1Point.theta)
        );
        
        // Check against all points in the second plate
        for (size_t j = 0; j < plate2.boundaries.size(); j += 3) {
            if (j + 2 >= plate2.boundaries.size()) break;
            
            Math::SphericalCoord p2Point{
                plate2.boundaries[j],
                plate2.boundaries[j+1],
                plate2.boundaries[j+2]
            };
            
            // Convert to cartesian
            Math::Vector3 p2CartPoint(
                p2Point.radius * sin(p2Point.theta) * cos(p2Point.phi),
                p2Point.radius * sin(p2Point.theta) * sin(p2Point.phi),
                p2Point.radius * cos(p2Point.theta)
            );
            
            // Check if points are close enough to be considered colliding
            float distance = (p1CartPoint - p2CartPoint).length();
            if (distance < collisionThreshold) {
                // This is a collision point - add the midpoint to our list
                Math::Vector3 collisionPoint = (p1CartPoint + p2CartPoint) * 0.5f;
                collisionPoints.push_back(collisionPoint);
            }
        }
    }
    
    // If we found any collision points, return true
    return !collisionPoints.empty();
}

void PlateTectonics::UpdateConvergentBoundaries(float deltaYears) {
    // Simple model: boundaries that have existed longer get higher
    for (auto& boundary : convergentBoundaries) {
        // Age the boundary
        boundary.height += boundary.heightChange * deltaYears;
    }
}

void PlateTectonics::LoadPlatesFromGenerator(const Tools::PlateGenerator& generator) {
    // Convert generator plates to tectonic plates
    plates = generator.convertToTectonicPlates();
    
    // Initialize mantle forces for all plates
    for (const auto& plate : plates) {
        mantleForces[plate.id] = 0.0f;
    }
    
    // Clear any existing convergent boundaries
    convergentBoundaries.clear();
}

std::vector<PlateTectonics::BoundaryPoint> PlateTectonics::GetConvergentBoundaries() const {
    return convergentBoundaries;
}

// Add after the existing CalculateMantleForces method
void PlateTectonics::CalculateImprovedMantleForces() {
    // Clear previous forces
    mantleForcesWithDirection.clear();
    
    // Create a more realistic mantle convection pattern
    // Based on Earth's major convection cells
    std::vector<MantelCell> convectionCells = CreateEarthLikeConvectionPattern();
    
    // For each plate, calculate forces from all convection cells
    for (const auto& plate : plates) {
        float totalForce = 0.0f;
        Math::Vector3 forceDirection(0, 0, 0);
        
        // Sample multiple points on the plate to get better force distribution
        std::vector<Math::Vector3> samplePoints = SamplePlatePoints(plate, 10);
        
        for (const auto& point : samplePoints) {
            // Calculate influence from each convection cell
            for (const auto& cell : convectionCells) {
                // Convert boundary points from spherical to cartesian for easier calculations
                Math::Vector3 cartPoint;
                
                // We need to have some point on the plate to calculate distance
                // This is a simplification - in a real implementation we would
                // properly sample points from plate.boundaries
                if (plate.boundaries.size() >= 3) {
                    Math::SphericalCoord spherePoint{
                        plate.boundaries[0],
                        plate.boundaries[1],
                        plate.boundaries[2]
                    };
                    
                    // Convert to cartesian
                    cartPoint = Math::Vector3(
                        spherePoint.radius * sin(spherePoint.theta) * cos(spherePoint.phi),
                        spherePoint.radius * sin(spherePoint.theta) * sin(spherePoint.phi),
                        spherePoint.radius * cos(spherePoint.theta)
                    );
                } else {
                    // Fallback if no boundaries
                    cartPoint = Math::Vector3(0, 0, 0);
                }
                
                float distance = (cartPoint - cell.center).length();
                float influence = CalculateCellInfluence(distance, cell);
                
                // Apply force in cell's flow direction
                forceDirection += cell.flowDirection * influence;
                totalForce += influence;
            }
        }
        
        // Average the force and normalize direction
        if (!samplePoints.empty()) {
            totalForce /= samplePoints.size();
            
            if (forceDirection.lengthSquared() > 0) {
                forceDirection = forceDirection.normalized();
            } else {
                forceDirection = Math::Vector3(0, 0, 1); // Default direction
            }
        }
        
        // Thickness affects resistance to movement
        float thicknessModifier = 100.0f / plate.thickness;
        
        // Density affects buoyancy
        float densityModifier = 3000.0f / plate.density;
        
        // Calculate final force with realistic magnitude
        float appliedForce = totalForce * thicknessModifier * densityModifier;
        
        // Store force with direction
        mantleForcesWithDirection[plate.id] = {appliedForce, forceDirection};
        
        // Also store in the original magnitude-only format for compatibility
        mantleForces[plate.id] = appliedForce;
    }
}

std::vector<PlateTectonics::MantelCell> PlateTectonics::CreateEarthLikeConvectionPattern() {
    std::vector<MantelCell> cells;
    
    // Earth has several major convection cells
    // We'll create a simplified model with 8 major cells
    
    // Create convection cells with specific patterns
    // Hadley cells (tropical)
    cells.push_back({
        Math::Vector3(0, 1, 0),         // Center near equator
        Math::Vector3(0, 0.2f, 1),      // Flow upward and slightly north
        1.2f                            // Stronger convection
    });
    
    cells.push_back({
        Math::Vector3(0, -1, 0),        // Center near equator (south)
        Math::Vector3(0, -0.2f, 1),     // Flow upward and slightly south
        1.2f                            // Stronger convection
    });
    
    // Ferrel cells (mid-latitude)
    cells.push_back({
        Math::Vector3(0, 0.5f, 0.7f),   // Northern mid-latitude
        Math::Vector3(0, -0.2f, -0.1f), // Flow mostly southward
        0.8f                            // Medium convection
    });
    
    cells.push_back({
        Math::Vector3(0, -0.5f, -0.7f), // Southern mid-latitude
        Math::Vector3(0, 0.2f, 0.1f),   // Flow mostly northward
        0.8f                            // Medium convection
    });
    
    // Polar cells
    cells.push_back({
        Math::Vector3(0, 0, 1),         // North pole
        Math::Vector3(0, -0.1f, -1),    // Flow downward and slightly south
        0.7f                            // Weaker convection
    });
    
    cells.push_back({
        Math::Vector3(0, 0, -1),        // South pole
        Math::Vector3(0, 0.1f, 1),      // Flow upward and slightly north
        0.7f                            // Weaker convection
    });
    
    // Add a few more cells for east-west variations
    cells.push_back({
        Math::Vector3(1, 0, 0),         // Eastern cell
        Math::Vector3(-1, 0.1f, 0),     // Flow westward
        0.6f                            // Weaker lateral convection
    });
    
    cells.push_back({
        Math::Vector3(-1, 0, 0),        // Western cell
        Math::Vector3(1, -0.1f, 0),     // Flow eastward
        0.6f                            // Weaker lateral convection
    });
    
    return cells;
}

std::vector<Math::Vector3> PlateTectonics::SamplePlatePoints(const Plate& plate, int sampleCount) {
    std::vector<Math::Vector3> points;
    
    // Early exit if no boundaries to sample from
    if (plate.boundaries.size() < 3) {
        return points;
    }
    
    // Sample points from plate boundaries
    // This is a simplification - we're just taking regular samples from available coordinates
    for (int i = 0; i < sampleCount; i++) {
        // Find an index that's safely within the boundaries array
        size_t index = (i * 3) % (plate.boundaries.size() - 2);
        
        if (index + 2 >= plate.boundaries.size()) {
            continue;
        }
        
        // Extract spherical coordinates
        Math::SphericalCoord spherePoint{
            plate.boundaries[index],
            plate.boundaries[index + 1],
            plate.boundaries[index + 2]
        };
        
        // Convert to cartesian
        Math::Vector3 cartPoint(
            spherePoint.radius * sin(spherePoint.theta) * cos(spherePoint.phi),
            spherePoint.radius * sin(spherePoint.theta) * sin(spherePoint.phi),
            spherePoint.radius * cos(spherePoint.theta)
        );
        
        points.push_back(cartPoint);
    }
    
    return points;
}

float PlateTectonics::CalculateCellInfluence(float distance, const MantelCell& cell) {
    // Calculate cell influence based on distance
    // Closer points are more influenced
    const float maxDistance = 2.0f; // Maximum effect distance
    
    if (distance >= maxDistance) {
        return 0.0f;
    }
    
    // Exponential falloff
    float normalizedDistance = distance / maxDistance;
    float influence = std::exp(-normalizedDistance * 3.0f) * cell.strength;
    
    return influence * 1e15f; // Scale to appropriate magnitude
}

// Add this method to PlateTectonics.cpp
void PlateTectonics::ProcessSubductionZones(float deltaSeconds) {
    std::vector<BoundaryPoint> newSubductionPoints;
    
    // Process each convergent boundary
    for (const auto& boundary : convergentBoundaries) {
        if (boundary.divergenceRate >= 0) continue; // Skip non-convergent boundaries
        
        // Get the plates involved
        int plate1Idx = -1, plate2Idx = -1;
        for (size_t i = 0; i < plates.size(); i++) {
            if (plates[i].id == boundary.plate1Id) plate1Idx = static_cast<int>(i);
            if (plates[i].id == boundary.plate2Id) plate2Idx = static_cast<int>(i);
        }
        
        if (plate1Idx < 0 || plate2Idx < 0) continue;
        
        // Determine which plate subducts (the denser one)
        bool plate1Subducts = plates[plate1Idx].density > plates[plate2Idx].density;
        int subductingIdx = plate1Subducts ? plate1Idx : plate2Idx;
        int overridingIdx = plate1Subducts ? plate2Idx : plate1Idx;
        
        // Calculate subduction speed based on convergence rate and plate properties
        float subductionSpeed = std::abs(boundary.divergenceRate) * 
                               (plates[subductingIdx].density / 3000.0f);  // Normalized for typical density
        
        // Modify thickness of subducting plate (thinning at leading edge)
        plates[subductingIdx].thickness *= (1.0f - 0.02f * deltaSeconds / (365.25f * 24.0f * 3600.0f));
        
        // Increases in density as plate subducts (metamorphism)
        plates[subductingIdx].density += 5.0f * deltaSeconds / (365.25f * 24.0f * 3600.0f);
        
        // Create volcanic arc in overriding plate
        // Distance from trench proportional to subduction angle (steeper = closer)
        float arcDistance = 200.0f * (1.0f - plates[subductingIdx].density / 3500.0f);
        
        // Calculate arc position
        Math::Vector3 toPlateCenter = (plates[overridingIdx].boundaries.size() >= 3) ?
            Math::Vector3(plates[overridingIdx].boundaries[0], 
                         plates[overridingIdx].boundaries[1], 
                         plates[overridingIdx].boundaries[2]) - boundary.position :
            Math::Vector3(0, 0, 1);
        
        toPlateCenter = toPlateCenter.normalized();
        
        // Volcanic arc position
        Math::Vector3 arcPosition = boundary.position + toPlateCenter * arcDistance;
        
        // Calculate mountain height based on subduction duration
        float mountainHeight = 500.0f + 5000.0f * (1.0f - std::exp(-boundary.age / 10000000.0f));
        
        // Calculate trench depth based on subducting plate age and density
        float trenchDepth = -2000.0f - 8000.0f * (plates[subductingIdx].density - 3000.0f) / 500.0f 
                          * (1.0f - std::exp(-plates[subductingIdx].age / 5000000.0f));
        
        // Add trench point
        BoundaryPoint trenchPoint = boundary;
        trenchPoint.height = trenchDepth;
        trenchPoint.heightChange = -subductionSpeed * 0.01f; // Trench deepens with time
        newSubductionPoints.push_back(trenchPoint);
        
        // Add volcanic arc point
        BoundaryPoint arcPoint;
        arcPoint.position = arcPosition;
        arcPoint.height = mountainHeight;
        arcPoint.heightChange = subductionSpeed * 0.02f; // Arc grows with time
        arcPoint.plate1Id = boundary.plate1Id;
        arcPoint.plate2Id = boundary.plate2Id;
        arcPoint.divergenceRate = boundary.divergenceRate;
        arcPoint.age = boundary.age;
        newSubductionPoints.push_back(arcPoint);
        
        // Modify plate velocities based on slab pull and resistance
        float slabPullForce = plates[subductingIdx].density * plates[subductingIdx].thickness * 9.8f * 0.01f;
        float resistiveForce = mantleViscosity * subductionSpeed * 1e-20f;
        
        // Apply forces to angular velocity
        plates[subductingIdx].angularVelocity += (slabPullForce - resistiveForce) * 
                                              deltaSeconds / (plates[subductingIdx].thickness * plates[subductingIdx].density);
        
        // Overriding plate gets pushed/pulled depending on coupling
        float couplingFactor = 0.2f; // 20% of slab force transfers to overriding plate
        plates[overridingIdx].angularVelocity += slabPullForce * couplingFactor * 
                                              deltaSeconds / (plates[overridingIdx].thickness * plates[overridingIdx].density);
    }
    
    // Update convergent boundaries with new points
    convergentBoundaries = newSubductionPoints;
}

// Add to PlateTectonics.cpp
void PlateTectonics::ProcessSpreadingRidges(float deltaSeconds) {
    std::vector<BoundaryPoint> newDivergentPoints;
    
    // Process each divergent boundary
    for (const auto& boundary : divergentBoundaries) {
        if (boundary.divergenceRate <= 0) continue; // Skip non-divergent boundaries
        
        // Get the plates involved
        int plate1Idx = -1, plate2Idx = -1;
        for (size_t i = 0; i < plates.size(); i++) {
            if (plates[i].id == boundary.plate1Id) plate1Idx = static_cast<int>(i);
            if (plates[i].id == boundary.plate2Id) plate2Idx = static_cast<int>(i);
        }
        
        if (plate1Idx < 0 || plate2Idx < 0) continue;
        
        // Calculate spreading rate based on boundary divergence rate
        float spreadingRate = boundary.divergenceRate;
        
        // Create new oceanic crust at the spreading center
        // Both plates get thinner near the boundary (conservation of mass)
        float crustFormationRate = spreadingRate * 0.01f * deltaSeconds / (365.25f * 24.0f * 3600.0f);
        
        // Ridge elevation based on thermal expansion of new crust
        float ridgeElevation = -2000.0f + 2500.0f * std::exp(-boundary.age / 2000000.0f);
        
        // Update boundary point
        BoundaryPoint updatedPoint = boundary;
        updatedPoint.height = ridgeElevation;
        updatedPoint.heightChange = -0.5f * spreadingRate * 0.001f; // Ridge subsides as it cools
        updatedPoint.age += deltaSeconds;
        newDivergentPoints.push_back(updatedPoint);
        
        // Create flanking bathymetry (ocean floor getting deeper with age)
        for (int side = -1; side <= 1; side += 2) {
            if (side == 0) continue;
            
            // Skip center point
            for (int distanceStep = 1; distanceStep <= 5; distanceStep++) {
                float distance = distanceStep * 50000.0f; // 50km steps
                
                // Compute point on ridge flank
                // This is simplified - real implementation would calculate actual geodesic offset
                Math::Vector3 flankDir = boundary.position.cross(Math::Vector3(0, 0, 1)).normalized();
                if (side < 0) flankDir = flankDir * -1.0f;
                
                Math::Vector3 flankPos = boundary.position + flankDir * (distance / 6371000.0f); // Normalized by Earth radius
                flankPos = flankPos.normalized() * 6371000.0f; // Back to Earth radius
                
                // Age of crust increases with distance from ridge
                float crustAge = distance / spreadingRate;
                
                // GDH1 model for ocean floor depth vs age
                float depth = -2500.0f - 350.0f * std::sqrt(crustAge / 1000000.0f);
                
                // Create a boundary point for this section of seafloor
                BoundaryPoint seafloorPoint;
                seafloorPoint.position = flankPos;
                seafloorPoint.height = depth;
                seafloorPoint.heightChange = -0.1f * spreadingRate * 0.001f; // Slower subsidence with age
                seafloorPoint.plate1Id = (side < 0) ? boundary.plate1Id : boundary.plate2Id;
                seafloorPoint.plate2Id = -1; // Not a plate boundary
                seafloorPoint.divergenceRate = 0;
                seafloorPoint.age = crustAge;
                
                newDivergentPoints.push_back(seafloorPoint);
            }
        }
        
        // Apply ridge push force to plates
        float ridgePushForce = 2.0f * mantleViscosity * spreadingRate * 1e-19f;
        
        // Ridge push acts in opposite directions for each plate
        plates[plate1Idx].angularVelocity += ridgePushForce * 
                                          deltaSeconds / (plates[plate1Idx].thickness * plates[plate1Idx].density);
        plates[plate2Idx].angularVelocity -= ridgePushForce * 
                                          deltaSeconds / (plates[plate2Idx].thickness * plates[plate2Idx].density);
                                          
        // Make oceanic plates get thinner with spreading (conservation of mass)
        if (plates[plate1Idx].isOceanic) {
            plates[plate1Idx].thickness *= (1.0f - 0.005f * crustFormationRate);
        }
        if (plates[plate2Idx].isOceanic) {
            plates[plate2Idx].thickness *= (1.0f - 0.005f * crustFormationRate);
        }
    }
    
    // Update divergent boundaries with new points
    divergentBoundaries = newDivergentPoints;
}

// Add to PlateTectonics.cpp
void PlateTectonics::ProcessContinentalCollisions(float deltaSeconds) {
    std::vector<BoundaryPoint> newCollisionPoints;
    
    // Process each convergent boundary between continental plates
    for (const auto& boundary : convergentBoundaries) {
        if (boundary.divergenceRate >= 0) continue; // Skip non-convergent boundaries
        
        // Get the plates involved
        int plate1Idx = -1, plate2Idx = -1;
        for (size_t i = 0; i < plates.size(); i++) {
            if (plates[i].id == boundary.plate1Id) plate1Idx = static_cast<int>(i);
            if (plates[i].id == boundary.plate2Id) plate2Idx = static_cast<int>(i);
        }
        
        if (plate1Idx < 0 || plate2Idx < 0) continue;
        
        // Skip if either plate is oceanic (handled by subduction)
        if (plates[plate1Idx].isOceanic || plates[plate2Idx].isOceanic) continue;
        
        // Calculate collision intensity from convergence rate
        float collisionRate = std::abs(boundary.divergenceRate);
        
        // Both continental plates thicken at collision zone
        float thickeningRate = collisionRate * 0.02f * deltaSeconds / (365.25f * 24.0f * 3600.0f);
        plates[plate1Idx].thickness *= (1.0f + thickeningRate * 0.5f);
        plates[plate2Idx].thickness *= (1.0f + thickeningRate * 0.5f);
        
        // Calculate mountain height based on:
        // 1. Collision duration (boundary age)
        // 2. Thickness of the colliding plates
        // 3. Collision rate
        // 4. Isostatic compensation (buoyancy)
        
        float combinedThickness = plates[plate1Idx].thickness + plates[plate2Idx].thickness;
        float collisionDuration = boundary.age;
        float buoyancyFactor = (3300.0f - plates[plate1Idx].density) / 600.0f + 
                              (3300.0f - plates[plate2Idx].density) / 600.0f;
        
        // S-curve growth function for mountains
        float maxHeight = 5000.0f + combinedThickness * 20.0f * buoyancyFactor;
        float growthRate = collisionRate * 0.2f;
        float halfTime = 10000000.0f; // 10 million years to reach half height
        float currentHeight = maxHeight / (1.0f + std::exp(-(collisionDuration - halfTime) * growthRate / halfTime));
        
        // Rate of height change
        float heightChange = maxHeight * growthRate * std::exp(-(collisionDuration - halfTime) * growthRate / halfTime) /
                           (1.0f + std::exp(-(collisionDuration - halfTime) * growthRate / halfTime)) / 
                           (1.0f + std::exp(-(collisionDuration - halfTime) * growthRate / halfTime));
        
        // Apply height changes to collision boundary
        BoundaryPoint updatedPoint = boundary;
        updatedPoint.height = currentHeight;
        updatedPoint.heightChange = heightChange;
        updatedPoint.age += deltaSeconds;
        newCollisionPoints.push_back(updatedPoint);
        
        // Create a mountain range with multiple peaks
        const int RANGE_WIDTH = 5;
        for (int offset = -RANGE_WIDTH; offset <= RANGE_WIDTH; offset++) {
            if (offset == 0) continue; // Center point already added
            
            // Calculate offset position along collision boundary
            // This is simplified - real implementation would follow actual boundary geometry
            Math::Vector3 offsetDir = boundary.position.cross(Math::Vector3(0, 1, 0)).normalized();
            float offsetDist = offset * 50000.0f; // 50km between mountain peaks
            
            Math::Vector3 peakPos = boundary.position + offsetDir * (offsetDist / 6371000.0f);
            peakPos = peakPos.normalized() * 6371000.0f;
            
            // Height varies along range - slightly lower than central peak
            float peakHeightFactor = 1.0f - 0.1f * std::abs(offset) / RANGE_WIDTH;
            float peakHeight = currentHeight * peakHeightFactor;
            
            // Add mountain peak
            BoundaryPoint peakPoint;
            peakPoint.position = peakPos;
            peakPoint.height = peakHeight;
            peakPoint.heightChange = heightChange * peakHeightFactor;
            peakPoint.plate1Id = boundary.plate1Id;
            peakPoint.plate2Id = boundary.plate2Id;
            peakPoint.divergenceRate = boundary.divergenceRate;
            peakPoint.age = boundary.age;
            
            newCollisionPoints.push_back(peakPoint);
        }
        
        // Continental collision slows down both plates (high resistance)
        float resistiveForce = mantleViscosity * collisionRate * 1e-18f;
        
        plates[plate1Idx].angularVelocity *= (1.0f - resistiveForce * 
                                          deltaSeconds / (plates[plate1Idx].thickness * plates[plate1Idx].density));
        plates[plate2Idx].angularVelocity *= (1.0f - resistiveForce * 
                                          deltaSeconds / (plates[plate2Idx].thickness * plates[plate2Idx].density));
    }
    
    // Add new collision points to the convergent boundaries list
    convergentBoundaries.insert(convergentBoundaries.end(), newCollisionPoints.begin(), newCollisionPoints.end());
}

void PlateTectonics::ClassifyPlateBoundaries() {
    // Clear previous boundary classifications
    divergentBoundaries.clear();
    transformBoundaries.clear();
    
    // Reclassify convergent boundaries
    std::vector<BoundaryPoint> newConvergentBoundaries;
    
    // Detect all collisions
    std::vector<PlateCollision> collisions = DetectAllCollisions();
    
    // For each collision, create appropriate boundary points
    for (const auto& collision : collisions) {
        // Ensure indices are valid
        if (collision.plate1Index < 0 || collision.plate1Index >= static_cast<int>(plates.size()) ||
            collision.plate2Index < 0 || collision.plate2Index >= static_cast<int>(plates.size())) {
            continue;
        }
        
        const Plate& plate1 = plates[collision.plate1Index];
        const Plate& plate2 = plates[collision.plate2Index];
        
        // Calculate relative velocity
        float relativeVelocity = plate1.angularVelocity - plate2.angularVelocity;
        
        // Process each collision point
        for (const auto& point : collision.points) {
            BoundaryPoint boundary;
            boundary.position = point;
            boundary.plate1Id = plate1.id;
            boundary.plate2Id = plate2.id;
            boundary.age = 0; // New boundary point
            
            // Classify based on relative motion
            if (std::abs(relativeVelocity) < 0.00001f) {
                // Transform boundary
                boundary.divergenceRate = 0;
                boundary.slipRate = 0.01f; // Small default value
                boundary.height = 0;
                boundary.heightChange = 0;
                transformBoundaries.push_back(boundary);
            }
            else if (relativeVelocity > 0) {
                // Convergent boundary
                boundary.divergenceRate = -relativeVelocity * 1000.0f * 31536000.0f; // Convert to mm/year
                boundary.slipRate = 0;
                
                // Initial height depends on plate types
                if (plate1.isOceanic && plate2.isOceanic) {
                    // Ocean-ocean: depends on age difference
                    float ageDiff = std::abs(plate1.age - plate2.age);
                    boundary.height = -4000.0f - 2000.0f * (ageDiff / 100000000.0f);
                }
                else if (!plate1.isOceanic && !plate2.isOceanic) {
                    // Continent-continent: initial low mountain
                    boundary.height = 500.0f;
                }
                else {
                    // Ocean-continent: initial trench
                    boundary.height = -3000.0f;
                }
                
                // Height change rate proportional to convergence
                boundary.heightChange = std::abs(relativeVelocity) * 100.0f;
                
                newConvergentBoundaries.push_back(boundary);
            }
            else {
                // Divergent boundary
                boundary.divergenceRate = -relativeVelocity * 1000.0f * 31536000.0f; // Convert to mm/year (positive for divergent)
                boundary.slipRate = 0;
                
                // New divergent boundaries start slightly below sea level
                boundary.height = -2500.0f;
                boundary.heightChange = 0;
                
                divergentBoundaries.push_back(boundary);
            }
        }
    }
    
    // Keep existing convergent boundaries and add newly detected ones
    for (auto& boundary : convergentBoundaries) {
        // Age existing boundaries
        boundary.age += 100000; // Approximately 100,000 years between simulation frames
    }
    
    // Add new boundaries to the list
    convergentBoundaries.insert(convergentBoundaries.end(), 
                              newConvergentBoundaries.begin(), 
                              newConvergentBoundaries.end());
}

// Add this to PlateTectonics.cpp
void PlateTectonics::UpdatePlateAngularMomentum(float deltaSeconds) {
    // Apply forces to plates and update velocities
    for (auto& plate : plates) {
        // Get force vector for this plate (with safety check)
        if (mantleForcesWithDirection.find(plate.id) != mantleForcesWithDirection.end()) {
            ForceVector forceVec = mantleForcesWithDirection[plate.id];
            
            // Apply mantle force to angular velocity
            float mantleForce = forceVec.magnitude;
            plate.angularVelocity += mantleForce * deltaSeconds / (plate.thickness * plate.density);
        }
        
        // Apply viscous drag from mantle
        float dragCoefficient = 1.0f / mantleViscosity * 1e19f; // Scale factor
        plate.angularVelocity *= (1.0f - dragCoefficient * deltaSeconds);
    }
}

// Add this to PlateTectonics.cpp too
void PlateTectonics::ProcessTransformFaults(float deltaSeconds) {
    // Simple implementation for transform faults
    for (auto& boundary : transformBoundaries) {
        // Transform faults mostly just slide past each other with minimal elevation changes
        boundary.height *= 0.99f; // Gradual erosion of any elevated features
        
        // Minor changes to position based on slip rate
        // In a real implementation, this would involve complex lateral movements
        
        // Age the boundary
        boundary.age += deltaSeconds;
    }
}
} // namespace Simulation
} // namespace AeonTerra