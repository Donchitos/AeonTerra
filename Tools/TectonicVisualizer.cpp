#include "TectonicVisualizer.h"
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>

namespace AeonTerra {
namespace Tools {

// Helper function to create directory if it doesn't exist
void EnsureDirectoryExists(const std::string& directory) {
    std::filesystem::path dir(directory);
    if (!std::filesystem::exists(dir)) {
        std::filesystem::create_directories(dir);
    }
}

void TectonicVisualizer::GeneratePlateTectonicAnimation(
    const Simulation::PlateTectonics& tectonics,
    float timeStart,
    float timeEnd,
    float timeStep,
    const std::string& outputDir) {
    
    // Ensure output directory exists
    EnsureDirectoryExists(outputDir);
    
    // Clone the tectonics object so we can simulate forward
    Simulation::PlateTectonics simulationCopy = tectonics; // Assumes copy constructor works
    
    // Number of frames to generate
    int frameCount = static_cast<int>((timeStart - timeEnd) / timeStep) + 1;
    
    // Generate each frame
    for (int frame = 0; frame < frameCount; frame++) {
        float currentTime = timeStart - (frame * timeStep);
        
        // Format frame number with leading zeros
        char frameNumStr[10];
        sprintf(frameNumStr, "%04d", frame);
        
        // Generate filename for this frame
        std::string filename = outputDir + "/plate_animation_" + frameNumStr + ".png";
        
        // Generate heightmap for current state
        HeightmapGenerator heightmapGen;
        HeightmapGenerator::Options options;
        options.resolution = 1024; // Higher resolution for animation
        
        std::vector<float> heightmap = heightmapGen.GenerateHeightmap(simulationCopy, options);
        
        // Save visualization to file
        heightmapGen.SaveShadedReliefMap(filename, heightmap, options.resolution);
        
        // Also save plate boundaries
        std::string plateFilename = outputDir + "/plate_bounds_" + frameNumStr + ".png";
        heightmapGen.SavePlateMap(plateFilename, simulationCopy, options.resolution);
        
        // Simulate forward to next frame
        simulationCopy.SimulateFrame(timeStep * 1000000.0f); // Convert to years
    }
}

void TectonicVisualizer::GenerateCrustMantelSection(
    const Simulation::PlateTectonics& tectonics,
    float latitude,
    float longitude,
    float depth,
    const std::string& outputFile) {
    
    // Create directory if needed
    std::string directory = outputFile.substr(0, outputFile.find_last_of("/\\"));
    EnsureDirectoryExists(directory);
    
    // Image dimensions
    int width = 1200;
    int height = 600;
    
    // Convert latitude and longitude to radians
    float latRad = latitude * M_PI / 180.0f;
    float lonRad = longitude * M_PI / 180.0f;
    
    // Calculate section start and end points
    const float planetRadius = 6371.0f; // km (Earth radius)
    
    // Open file for binary writing (PPM format)
    std::ofstream file(outputFile, std::ios::binary);
    if (!file.is_open()) {
        return;
    }
    
    // PPM header (P6 format)
    file << "P6\n" << width << " " << height << "\n255\n";
    
    // Get plate data
    const auto& plates = tectonics.GetPlates();
    const auto& mantleCells = tectonics.GetMantleCells();
    
    // Create the cross-section image
    for (int y = 0; y < height; y++) {
        // Convert y to depth in planet (0 at surface, height-1 at maximum depth)
        float depthFraction = static_cast<float>(y) / (height - 1);
        float currentDepth = depth * depthFraction;
        float radius = planetRadius - currentDepth;
        
        for (int x = 0; x < width; x++) {
            // Determine color based on depth/layer
            unsigned char r, g, b;
            
            if (radius >= planetRadius - 35.0f) {
                // Crust layer
                r = 150;
                g = 120;
                b = 90;
            }
            else if (radius >= planetRadius - 100.0f) {
                // Upper mantle
                r = 180;
                g = 160;
                b = 50;
            }
            else if (radius >= planetRadius - 660.0f) {
                // Lower mantle
                r = 200;
                g = 100;
                b = 0;
            }
            else if (radius >= planetRadius - 2900.0f) {
                // Deep mantle
                r = 160;
                g = 0;
                b = 0;
            }
            else if (radius >= planetRadius - 5100.0f) {
                // Outer core
                r = 230;
                g = 150;
                b = 0;
            }
            else {
                // Inner core
                r = 255;
                g = 220;
                b = 0;
            }
            
            // Write RGB pixel
            file.put(static_cast<char>(r));
            file.put(static_cast<char>(g));
            file.put(static_cast<char>(b));
        }
    }
}

void TectonicVisualizer::GenerateStressMap(
    const Simulation::PlateTectonics& tectonics,
    const std::string& outputFile) {
    
    // Create directory if needed
    std::string directory = outputFile.substr(0, outputFile.find_last_of("/\\"));
    EnsureDirectoryExists(directory);
    
    // Image dimensions 
    int width = 2048;
    int height = 1024;
    
    // Open file for binary writing (PPM format)
    std::ofstream file(outputFile, std::ios::binary);
    if (!file.is_open()) {
        return;
    }
    
    // PPM header (P6 format)
    file << "P6\n" << width << " " << height << "\n255\n";
    
    // Get tectonic data
    const auto& plates = tectonics.GetPlates();
    const auto& boundaries = tectonics.GetConvergentBoundaries();
    
    // Create temporary stress map array
    std::vector<float> stressMap(width * height, 0.0f);
    
    // Compute stress field from plate boundaries
    // Higher stress near boundaries, especially convergent ones
    for (const auto& boundary : boundaries) {
        // Convert 3D position to spherical, then to image coordinates
        Math::Vector3 pos = boundary.position;
        float len = pos.length();
        if (len < 0.001f) continue;
        
        // Convert to spherical coordinates
        float theta = std::acos(pos.z / len);
        float phi = std::atan2(pos.y, pos.x);
        if (phi < 0) phi += 2 * M_PI;
        
        // Convert to image coordinates
        int x = static_cast<int>((phi / (2 * M_PI)) * width) % width;
        int y = static_cast<int>((theta / M_PI) * height) % height;
        
        // Calculate stress based on boundary type
        float stress = 5.0f; // Base stress
        
        // Add stress to nearby pixels with distance falloff
        const int radius = 50; // Influence radius
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {
                int nx = (x + dx + width) % width;
                int ny = (y + dy + height) % height;
                
                // Distance-based falloff
                float dist = std::sqrt(dx*dx + dy*dy);
                if (dist > radius) continue;
                
                float influence = (radius - dist) / radius;
                stressMap[ny * width + nx] += stress * influence;
            }
        }
    }
    
    // Find max stress for normalization
    float maxStress = *std::max_element(stressMap.begin(), stressMap.end());
    if (maxStress < 0.001f) maxStress = 0.001f; // Avoid division by zero
    
    // Create the stress map image
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = y * width + x;
            float stress = stressMap[index] / maxStress;
            
            // Convert to color (blue-white-red gradient)
            unsigned char r, g, b;
            
            if (stress < 0.5f) {
                // Blue to white
                float t = stress * 2.0f;
                r = static_cast<unsigned char>(t * 255);
                g = static_cast<unsigned char>(t * 255);
                b = 255;
            } else {
                // White to red
                float t = (stress - 0.5f) * 2.0f;
                r = 255;
                g = static_cast<unsigned char>((1.0f - t) * 255);
                b = static_cast<unsigned char>((1.0f - t) * 255);
            }
            
            // Write RGB pixel
            file.put(static_cast<char>(r));
            file.put(static_cast<char>(g));
            file.put(static_cast<char>(b));
        }
    }
}

} // namespace Tools
} // namespace AeonTerra