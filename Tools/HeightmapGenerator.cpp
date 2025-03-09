#include "HeightmapGenerator.h"
#include "../Core/Mathematics/GeoMath.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <limits>

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Make sure we're building with correct export macros
#ifdef AEONTERRA_CORE_STATIC
#define AEONTERRA_API
#endif

namespace AeonTerra {
namespace Tools {

HeightmapGenerator::HeightmapGenerator() {
    // No initialization needed at this time
}

std::vector<float> HeightmapGenerator::GenerateHeightmap(
    const Simulation::PlateTectonics& tectonics,
    const Options& options) {
    
    // Create empty heightmap
    int size = options.resolution * options.resolution;
    std::vector<float> heightmap(size, 0.0f);
    
    // Get plate data
    const auto& plates = tectonics.GetPlates();
    
    // Generate heightmap points
    #pragma omp parallel for
    for (int y = 0; y < options.resolution; ++y) {
        for (int x = 0; x < options.resolution; ++x) {
            // Convert from grid coordinates to spherical coordinates
            float u = static_cast<float>(x) / (options.resolution - 1);
            float v = static_cast<float>(y) / (options.resolution - 1);
            
            // Map to spherical coordinates
            float theta = v * M_PI;          // 0 to PI (north to south)
            float phi = u * 2.0f * M_PI;     // 0 to 2PI (longitude)
            
            // Convert to 3D point on unit sphere
            Math::Vector3 point(
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)
            );
            
            // Calculate height at this point
            float height = CalculateHeight(point, tectonics, options);
            
            // Store in heightmap
            int index = y * options.resolution + x;
            heightmap[index] = height;
        }
    }
    
    return heightmap;
}

float HeightmapGenerator::CalculateHeight(
    const Math::Vector3& point,
    const Simulation::PlateTectonics& tectonics,
    const Options& options) {
    
    const auto& plates = tectonics.GetPlates();
    
    // Find the plate this point belongs to
    // In a real implementation, we would use spatial partitioning for efficiency
    int closestPlateIndex = -1;
    float minDistance = std::numeric_limits<float>::max();
    
    for (size_t i = 0; i < plates.size(); ++i) {
        // Calculate distance to plate center
        // Note: This is a simplified approach. In a complete implementation,
        // we would check point-in-polygon for the plate boundaries
        
        // For this simplified version, we'll just use the first three values
        // in the boundaries vector as the plate center (if available)
        if (plates[i].boundaries.size() >= 3) {
            Math::SphericalCoord plateCenter{
                plates[i].boundaries[0],
                plates[i].boundaries[1],
                plates[i].boundaries[2]
            };
            
            Math::SphericalCoord pointSpherical{
                1.0f, // Unit radius
                std::acos(point.z),
                std::atan2(point.y, point.x)
            };
            
            float distance = Math::GeoMath::GreatCircleDistance(plateCenter, pointSpherical);
            
            if (distance < minDistance) {
                minDistance = distance;
                closestPlateIndex = i;
            }
        }
    }
    
    // Base height depends on plate density
    // Lower density plates form continents, higher density plates form ocean basins
    float baseHeight = 0.0f;
    
    if (closestPlateIndex >= 0) {
        const auto& plate = plates[closestPlateIndex];
        
        // Higher density = deeper ocean
        // Lower density = higher continent
        float densityFactor = (3300.0f - plate.density) / 600.0f; // Normalized to roughly -0.5 to 0.5
        
        if (densityFactor > 0) {
            // Continental plate
            baseHeight = options.continentHeight * densityFactor;
        } else {
            // Oceanic plate
            baseHeight = options.oceanDepth * -densityFactor;
        }
        
        // Thicker plates have higher elevation
        float thicknessFactor = (plate.thickness - 50.0f) / 50.0f; // Normalized to roughly -0.5 to 0.5
        baseHeight += thicknessFactor * 500.0f;
        
        // Add mountain ranges near plate boundaries
        float boundaryDistance = DistanceToPlateBoundary(point, tectonics);
        float boundaryFactor = std::exp(-boundaryDistance * 10.0f);
        
        // Only add mountains on continental plates
        if (densityFactor > 0) {
            baseHeight += boundaryFactor * options.mountainScale;
        } else {
            // For oceanic plates, create trenches
            baseHeight -= boundaryFactor * options.mountainScale * 0.5f;
        }
    }
    
    // Add some noise for natural-looking terrain
    float noise = GenerateNoise(point.x * 10.0f, point.y * 10.0f, options.noiseScale);
    baseHeight += noise;
    
    return baseHeight;
}

float HeightmapGenerator::DistanceToPlateBoundary(
    const Math::Vector3& point,
    const Simulation::PlateTectonics& tectonics) {
    
    // Simplified method - in a full implementation this would use
    // proper polygon edge distance calculations
    return 0.1f; // Placeholder value
}

float HeightmapGenerator::GenerateNoise(float x, float y, float scale) {
    // Simple Perlin-like noise function
    // In a real implementation, we would use a proper noise library
    
    x = std::fmod(x, 1.0f) * 6.28f;
    y = std::fmod(y, 1.0f) * 6.28f;
    
    return (std::sin(x) * std::cos(y) * 0.5f + 0.5f) * scale;
}

bool HeightmapGenerator::SaveHeightmap(
    const std::string& filename,
    const std::vector<float>& heightmap,
    int resolution) {
    
    // Open file for binary writing
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    // Write raw heightmap data
    file.write(reinterpret_cast<const char*>(heightmap.data()), 
               heightmap.size() * sizeof(float));
    
    return true;
}

} // namespace Tools
} // namespace AeonTerra