#include "HeightmapGenerator.h"
#include "../Core/Mathematics/GeoMath.h"
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>
#include <limits>
#include <tuple>

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

HeightmapGenerator::HeightmapGenerator() : m_seed(12345) {
    // Initialize with default seed
}

void HeightmapGenerator::SetSeed(uint32_t seed) {
    m_seed = seed;
}

std::vector<float> HeightmapGenerator::GenerateHeightmap(
    const Simulation::PlateTectonics& tectonics,
    const Options& options) {
    
    // Create empty heightmap
    int size = options.resolution * options.resolution;
    std::vector<float> heightmap(size, 0.0f);
    
    // Get plate data
    const auto& plates = tectonics.GetPlates();
    
    // Get convergent boundaries for mountain/trench visualization
    const auto& boundaries = tectonics.GetConvergentBoundaries();
    
    // Create noise generator with our seed
    std::mt19937 noiseGen(m_seed);
    std::uniform_real_distribution<float> noiseDist(-1.0f, 1.0f);
    
    // Generate fractal noise lookup grid for terrain detail
    const int noiseGridSize = 64;
    std::vector<float> noiseGrid(noiseGridSize * noiseGridSize);
    for (int i = 0; i < static_cast<int>(noiseGrid.size()); i++) {
        noiseGrid[i] = noiseDist(noiseGen);
    }
    
    // Generate heightmap points
    #pragma omp parallel for
    for (int y = 0; y < options.resolution; ++y) {
        for (int x = 0; x < options.resolution; ++x) {
            // Convert from grid coordinates to spherical coordinates
            float u = static_cast<float>(x) / (options.resolution - 1);
            float v = static_cast<float>(y) / (options.resolution - 1);
            
            // Map to spherical coordinates
            float theta = v * static_cast<float>(M_PI);          // 0 to PI (north to south)
            float phi = u * 2.0f * static_cast<float>(M_PI);     // 0 to 2PI (longitude)
            
            // Convert to 3D point on unit sphere
            Math::Vector3 point(
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)
            );
            
            // Calculate base height at this point based on plate characteristics
            float height = CalculateBaseHeight(point, tectonics, options);
            
            // Add height from convergent boundaries
            height += CalculateBoundaryHeight(point, boundaries, options);
            
            // Add fractal noise for detailed terrain
            height += CalculateFractalNoise(point, noiseGrid, noiseGridSize, options.noiseScale);
            
            // Store in heightmap
            int index = y * options.resolution + x;
            heightmap[index] = height;
        }
    }
    
    // Normalize heightmap to ensure it has a good range
    NormalizeHeightmap(heightmap, options.oceanDepth, options.continentHeight + options.mountainScale);
    
    return heightmap;
}

float HeightmapGenerator::CalculateBaseHeight(
    const Math::Vector3& point,
    const Simulation::PlateTectonics& tectonics,
    const Options& options) {
    
    const auto& plates = tectonics.GetPlates();
    
    // Find the plate this point belongs to
    int closestPlateIndex = -1;
    float minDistance = std::numeric_limits<float>::max();
    
    for (size_t i = 0; i < plates.size(); ++i) {
        // Calculate distance to plate center
        // Note: This is a simplified approach.
        
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
                closestPlateIndex = static_cast<int>(i);
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
    }
    
    return baseHeight;
}

float HeightmapGenerator::CalculateBoundaryHeight(
    const Math::Vector3& point,
    const std::vector<Simulation::PlateTectonics::BoundaryPoint>& boundaries,
    const Options& options) {
    
    float totalBoundaryHeight = 0.0f;
    
    // Check against all boundary points
    for (const auto& boundary : boundaries) {
        // Calculate distance to boundary point
        float distance = (point - boundary.position).length();
        
        // Use exponential falloff to create mountain ranges and trenches
        const float falloffFactor = 5.0f; // Controls width of mountain range
        float influenceFactor = std::exp(-distance * falloffFactor);
        
        // Apply the boundary height with falloff
        totalBoundaryHeight += boundary.height * influenceFactor;
    }
    
    return totalBoundaryHeight;
}

float HeightmapGenerator::CalculateFractalNoise(
    const Math::Vector3& point, 
    const std::vector<float>& noiseGrid, 
    int gridSize,
    float scale) {
    
    // Generate fractal noise by sampling multiple octaves
    float noise = 0.0f;
    float amplitude = 1.0f;
    float frequency = 1.0f;
    float totalAmplitude = 0.0f;
    
    const int octaves = 4;
    const float persistence = 0.5f;
    
    for (int i = 0; i < octaves; i++) {
        // Sample noise by using the point's coordinates
        float nx = std::fmod(std::abs(point.x * frequency * 10.0f), 1.0f) * static_cast<float>(gridSize - 1);
        float ny = std::fmod(std::abs(point.y * frequency * 10.0f), 1.0f) * static_cast<float>(gridSize - 1);
        
        // Bilinear interpolation for smoother noise
        int x0 = static_cast<int>(nx);
        int y0 = static_cast<int>(ny);
        int x1 = (x0 + 1) % gridSize;
        int y1 = (y0 + 1) % gridSize;
        
        float tx = nx - static_cast<float>(x0);
        float ty = ny - static_cast<float>(y0);
        
        float c00 = noiseGrid[y0 * gridSize + x0];
        float c10 = noiseGrid[y0 * gridSize + x1];
        float c01 = noiseGrid[y1 * gridSize + x0];
        float c11 = noiseGrid[y1 * gridSize + x1];
        
        // Bilinear interpolation
        float nx0 = c00 * (1.0f - tx) + c10 * tx;
        float nx1 = c01 * (1.0f - tx) + c11 * tx;
        float n = nx0 * (1.0f - ty) + nx1 * ty;
        
        // Add to total noise
        noise += n * amplitude;
        
        // Update for next octave
        totalAmplitude += amplitude;
        amplitude *= persistence;
        frequency *= 2.0f;
    }
    
    // Normalize and scale
    noise /= totalAmplitude;
    return noise * scale;
}

void HeightmapGenerator::NormalizeHeightmap(
    std::vector<float>& heightmap,
    float minValue,
    float maxValue) {
    
    // Find current min and max
    float currentMin = *std::min_element(heightmap.begin(), heightmap.end());
    float currentMax = *std::max_element(heightmap.begin(), heightmap.end());
    
    // Avoid division by zero
    if (std::abs(currentMax - currentMin) < 0.0001f) {
        return;
    }
    
    // Normalize to desired range
    float scale = (maxValue - minValue) / (currentMax - currentMin);
    
    for (auto& height : heightmap) {
        height = minValue + (height - currentMin) * scale;
    }
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

bool HeightmapGenerator::SaveHeightmapPNG(
    const std::string& filename,
    const std::vector<float>& heightmap,
    int resolution) {
    
    // Open file for binary writing
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    // Writing a simple PPM format image (P6 format)
    // Header
    file << "P6\n" << resolution << " " << resolution << "\n255\n";
    
    // Convert heightmap values to RGB colors with enhanced visibility
    for (size_t i = 0; i < heightmap.size(); ++i) {
        float height = heightmap[i];
        
        // Determine color based on height with more distinctive color scheme
        unsigned char r, g, b;
        
        if (height < -6000) { // Deep ocean
            r = 0;
            g = 20;
            b = 120;
        } 
        else if (height < -4000) { // Ocean
            r = 0;
            g = 50;
            b = 170;
        }
        else if (height < -2000) { // Ocean
            r = 0;
            g = 90;
            b = 220;
        }
        else if (height < 0) { // Shallow water
            r = 10;
            g = static_cast<unsigned char>(120 + (height + 2000) / 2000.0f * 135.0f);
            b = 250;
        }
        else if (height < 500) { // Plains
            r = static_cast<unsigned char>(50 + height / 500.0f * 130.0f);
            g = static_cast<unsigned char>(180 + height / 500.0f * 70.0f);
            b = static_cast<unsigned char>(50);
        }
        else if (height < 1500) { // Hills
            r = 180;
            g = 160;
            b = 50;
        }
        else if (height < 3000) { // Mountains
            r = 180;
            g = 150;
            b = 120;
        }
        else if (height < 5000) { // High Mountains
            r = 200;
            g = 170;
            b = 170;
        }
        else { // Highest Peaks
            r = 255;
            g = 250;
            b = 250;
        }
        
        // Add contour lines for emphasis
        float contourInterval = 500.0f; // Interval between contour lines
        float contourWidth = 100.0f; // Width of contour lines
        float moduloHeight = std::fmod(std::abs(height), contourInterval);
        
        if (moduloHeight < contourWidth || moduloHeight > contourInterval - contourWidth) {
            // Darken color to create contour line effect
            r = static_cast<unsigned char>(r * 0.7f);
            g = static_cast<unsigned char>(g * 0.7f);
            b = static_cast<unsigned char>(b * 0.7f);
        }
        
        // Write the RGB triple
        file.put(static_cast<char>(r));
        file.put(static_cast<char>(g));
        file.put(static_cast<char>(b));
    }
    
    return true;
}

bool HeightmapGenerator::SaveComparisonImage(
    const std::string& filename,
    const std::vector<float>& beforeHeightmap,
    const std::vector<float>& afterHeightmap,
    int resolution) {
    
    // Open file for binary writing
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    // We'll create a side-by-side image, so the width is doubled
    file << "P6\n" << (resolution * 2) << " " << resolution << "\n255\n";
    
    // Function to convert height to RGB
    auto heightToRGB = [](float height) -> std::tuple<unsigned char, unsigned char, unsigned char> {
        unsigned char r, g, b;
        
        if (height < -6000) { // Deep ocean
            r = 0;
            g = 20;
            b = 120;
        } 
        else if (height < -4000) { // Ocean
            r = 0;
            g = 50;
            b = 170;
        }
        else if (height < -2000) { // Ocean
            r = 0;
            g = 90;
            b = 220;
        }
        else if (height < 0) { // Shallow water
            r = 10;
            g = static_cast<unsigned char>(120 + (height + 2000) / 2000.0f * 135.0f);
            b = 250;
        }
        else if (height < 500) { // Plains
            r = static_cast<unsigned char>(50 + height / 500.0f * 130.0f);
            g = static_cast<unsigned char>(180 + height / 500.0f * 70.0f);
            b = static_cast<unsigned char>(50);
        }
        else if (height < 1500) { // Hills
            r = 180;
            g = 160;
            b = 50;
        }
        else if (height < 3000) { // Mountains
            r = 180;
            g = 150;
            b = 120;
        }
        else if (height < 5000) { // High Mountains
            r = 200;
            g = 170;
            b = 170;
        }
        else { // Highest Peaks
            r = 255;
            g = 250;
            b = 250;
        }
        
        return std::make_tuple(r, g, b);
    };
    
    // Generate side-by-side comparison
    for (int y = 0; y < resolution; ++y) {
        for (int x = 0; x < resolution; ++x) {
            int index = y * resolution + x;
            
            // Draw before image on the left
            if (index < static_cast<int>(beforeHeightmap.size())) {
                auto [r, g, b] = heightToRGB(beforeHeightmap[index]);
                file.put(static_cast<char>(r));
                file.put(static_cast<char>(g));
                file.put(static_cast<char>(b));
            }
            else {
                // Default black if out of bounds
                file.put(static_cast<char>(0));
                file.put(static_cast<char>(0));
                file.put(static_cast<char>(0));
            }
        }
        
        // Draw after image on the right
        for (int x = 0; x < resolution; ++x) {
            int index = y * resolution + x;
            
            if (index < static_cast<int>(afterHeightmap.size())) {
                auto [r, g, b] = heightToRGB(afterHeightmap[index]);
                file.put(static_cast<char>(r));
                file.put(static_cast<char>(g));
                file.put(static_cast<char>(b));
            }
            else {
                // Default black if out of bounds
                file.put(static_cast<char>(0));
                file.put(static_cast<char>(0));
                file.put(static_cast<char>(0));
            }
        }
    }
    
    return true;
}

bool HeightmapGenerator::SaveShadedReliefMap(
    const std::string& filename,
    const std::vector<float>& heightmap,
    int resolution) {
    
    // Open file for binary writing
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    // Writing a simple PPM format image (P6 format)
    // Header
    file << "P6\n" << resolution << " " << resolution << "\n255\n";
    
    // First pass: create a normal map from the heightmap
    std::vector<Math::Vector3> normals(heightmap.size(), Math::Vector3(0, 0, 1));
    
    // Calculate normals using finite differences
    for (int y = 1; y < resolution - 1; ++y) {
        for (int x = 1; x < resolution - 1; ++x) {
            int index = y * resolution + x;
            
            // Sample neighbors
            float h = heightmap[index]; // Center
            float hL = heightmap[index - 1]; // Left
            float hR = heightmap[index + 1]; // Right
            float hT = heightmap[(y - 1) * resolution + x]; // Top
            float hB = heightmap[(y + 1) * resolution + x]; // Bottom
            
            // Calculate partial derivatives
            float dhdx = (hR - hL) / 2.0f;
            float dhdy = (hB - hT) / 2.0f;
            
            // Surface normal is cross product of tangent vectors
            Math::Vector3 normal(
                -dhdx,
                -dhdy,
                1.0f
            );
            
            // Normalize
            normals[index] = normal.normalized();
        }
    }
    
    // Second pass: apply lighting to create a shaded relief map
    // Light direction (from upper left)
    Math::Vector3 lightDir(-1.0f, -1.0f, 1.0f);
    lightDir = lightDir.normalized();
    
    for (size_t i = 0; i < heightmap.size(); ++i) {
        float height = heightmap[i];
        
        // Determine base color based on height
        unsigned char r, g, b;
        
        if (height < -6000) { // Deep ocean
            r = 0;
            g = 20;
            b = 120;
        } 
        else if (height < -4000) { // Ocean
            r = 0;
            g = 50;
            b = 170;
        }
        else if (height < -2000) { // Ocean
            r = 0;
            g = 90;
            b = 220;
        }
        else if (height < 0) { // Shallow water
            r = 10;
            g = static_cast<unsigned char>(120 + (height + 2000) / 2000.0f * 135.0f);
            b = 250;
        }
        else if (height < 500) { // Plains
            r = static_cast<unsigned char>(50 + height / 500.0f * 130.0f);
            g = static_cast<unsigned char>(180 + height / 500.0f * 70.0f);
            b = static_cast<unsigned char>(50);
        }
        else if (height < 1500) { // Hills
            r = 180;
            g = 160;
            b = 50;
        }
        else if (height < 3000) { // Mountains
            r = 180;
            g = 150;
            b = 120;
        }
        else if (height < 5000) { // High Mountains
            r = 200;
            g = 170;
            b = 170;
        }
        else { // Highest Peaks
            r = 255;
            g = 250;
            b = 250;
        }
        
        // Apply lighting
        float diffuse = std::max(0.0f, normals[i].dot(lightDir));
        float ambient = 0.3f;
        float lighting = ambient + diffuse * (1.0f - ambient);
        
        // Apply lighting to color
        r = static_cast<unsigned char>(std::min(255.0f, r * lighting));
        g = static_cast<unsigned char>(std::min(255.0f, g * lighting));
        b = static_cast<unsigned char>(std::min(255.0f, b * lighting));
        
        // Write the RGB triple
        file.put(static_cast<char>(r));
        file.put(static_cast<char>(g));
        file.put(static_cast<char>(b));
    }
    
    return true;
}

bool HeightmapGenerator::SavePlateMap(
    const std::string& filename,
    const Simulation::PlateTectonics& tectonics,
    int resolution) {
    
    // Open file for binary writing
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    // Get plates and boundaries
    const auto& plates = tectonics.GetPlates();
    const auto& boundaries = tectonics.GetConvergentBoundaries();
    
    // Create a color mapping for plates (up to 20 distinct colors)
    const int NUM_COLORS = 20;
    struct RGBColor { unsigned char r, g, b; };
    std::vector<RGBColor> plateColors(NUM_COLORS);
    
    // Generate distinct colors for plates
    for (int i = 0; i < NUM_COLORS; i++) {
        // HSV to RGB conversion for evenly spaced hues
        float h = static_cast<float>(i) / NUM_COLORS;
        float s = 0.7f;
        float v = 0.9f;
        
        float c = v * s;
        float x = c * (1.0f - std::abs(std::fmod(h * 6.0f, 2.0f) - 1.0f));
        float m = v - c;
        
        float r, g, b;
        if (h < 1.0f/6.0f) { r = c; g = x; b = 0; }
        else if (h < 2.0f/6.0f) { r = x; g = c; b = 0; }
        else if (h < 3.0f/6.0f) { r = 0; g = c; b = x; }
        else if (h < 4.0f/6.0f) { r = 0; g = x; b = c; }
        else if (h < 5.0f/6.0f) { r = x; g = 0; b = c; }
        else { r = c; g = 0; b = x; }
        
        plateColors[i] = {
            static_cast<unsigned char>((r + m) * 255),
            static_cast<unsigned char>((g + m) * 255),
            static_cast<unsigned char>((b + m) * 255)
        };
    }
    
    // Writing a simple PPM format image (P6 format)
    // Header
    file << "P6\n" << resolution << " " << resolution << "\n255\n";
    
    // Generate plate map
    std::vector<int> plateMap(resolution * resolution, -1);
    
    // Assign each point to nearest plate
    for (int y = 0; y < resolution; ++y) {
        for (int x = 0; x < resolution; ++x) {
            // Convert from grid coordinates to spherical coordinates
            float u = static_cast<float>(x) / (resolution - 1);
            float v = static_cast<float>(y) / (resolution - 1);
            
            // Map to spherical coordinates
            float theta = v * static_cast<float>(M_PI);          // 0 to PI (north to south)
            float phi = u * 2.0f * static_cast<float>(M_PI);     // 0 to 2PI (longitude)
            
            // Convert to 3D point on unit sphere
            Math::Vector3 point(
                std::sin(theta) * std::cos(phi),
                std::sin(theta) * std::sin(phi),
                std::cos(theta)
            );
            
            // Find closest plate
            int closestPlate = -1;
            float minDistance = std::numeric_limits<float>::max();
            
            for (size_t i = 0; i < plates.size(); ++i) {
                // Check if plate has boundaries
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
                        closestPlate = plates[i].id % NUM_COLORS;
                    }
                }
            }
            
            // Store plate ID
            int index = y * resolution + x;
            plateMap[index] = closestPlate;
        }
    }
    
    // Draw plate map with boundaries highlighted
    for (int y = 0; y < resolution; ++y) {
        for (int x = 0; x < resolution; ++x) {
            int index = y * resolution + x;
            int plateId = plateMap[index];
            
            // Default color (black for unassigned)
            unsigned char r = 0, g = 0, b = 0;
            
            // If assigned to a plate, use the plate color
            if (plateId >= 0 && plateId < NUM_COLORS) {
                r = plateColors[plateId].r;
                g = plateColors[plateId].g;
                b = plateColors[plateId].b;
            }
            
            // Check if this is a boundary by looking at neighbors
            bool isBoundary = false;
            
            if (x > 0 && plateMap[index - 1] != plateId) isBoundary = true;
            if (x < resolution - 1 && plateMap[index + 1] != plateId) isBoundary = true;
            if (y > 0 && plateMap[index - resolution] != plateId) isBoundary = true;
            if (y < resolution - 1 && plateMap[index + resolution] != plateId) isBoundary = true;
            
            // Draw boundaries in black
            if (isBoundary) {
                r = 0;
                g = 0;
                b = 0;
            }
            
            // Write the RGB triple
            file.put(static_cast<char>(r));
            file.put(static_cast<char>(g));
            file.put(static_cast<char>(b));
        }
    }
    
    return true;
}

} // namespace Tools
} // namespace AeonTerra