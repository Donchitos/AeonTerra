#include "Core/Simulation/PlateTectonics.h"
#include "Tools/WorldBuilder/PlateGenerator.h"
#include "Tools/HeightmapGenerator.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>

// Simple PPM image writer for visualization
void writePPM(const std::string& filename, const std::vector<float>& heightmap, int width, int height) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    
    // PPM header
    file << "P3\n" << width << " " << height << "\n255\n";
    
    // Find min/max for normalization
    float minHeight = *std::min_element(heightmap.begin(), heightmap.end());
    float maxHeight = *std::max_element(heightmap.begin(), heightmap.end());
    float range = maxHeight - minHeight;
    
    // Write pixel data
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int index = y * width + x;
            float normalizedHeight = (heightmap[index] - minHeight) / range;
            
            // Color based on height
            int r, g, b;
            
            if (heightmap[index] < 0) {
                // Ocean - deep blue to light blue
                float oceanDepth = -heightmap[index] / 8000.0f;
                oceanDepth = std::min(oceanDepth, 1.0f);
                r = static_cast<int>((1.0f - oceanDepth) * 100);
                g = static_cast<int>((1.0f - oceanDepth) * 150);
                b = static_cast<int>(128 + oceanDepth * 127);
            } else {
                // Land - green to brown to white
                float landHeight = heightmap[index] / 8000.0f;
                landHeight = std::min(landHeight, 1.0f);
                
                if (landHeight < 0.2f) {
                    // Lowlands (green)
                    r = static_cast<int>(100 + landHeight * 500);
                    g = static_cast<int>(200 - landHeight * 100);
                    b = static_cast<int>(100 + landHeight * 200);
                } else if (landHeight < 0.7f) {
                    // Mountains (brown)
                    float t = (landHeight - 0.2f) / 0.5f;
                    r = static_cast<int>(200 + t * 55);
                    g = static_cast<int>(180 - t * 100);
                    b = static_cast<int>(140 - t * 100);
                } else {
                    // Peaks (white)
                    float t = (landHeight - 0.7f) / 0.3f;
                    r = static_cast<int>(255);
                    g = static_cast<int>(255);
                    b = static_cast<int>(255);
                }
            }
            
            file << r << " " << g << " " << b << " ";
        }
        file << "\n";
    }
    
    file.close();
}

int main(int argc, char** argv) {
    std::cout << "AeonTerra Planet Generator Test" << std::endl;
    std::cout << "===============================" << std::endl;
    
    // Parameters
    int plateCount = 15;
    float planetRadius = 6371.0f; // Earth radius in km
    int resolution = 512;
    float simulationYears = 10.0f; // Million years
    int frameCount = 10;
    
    // Initialize random seed
    uint32_t seed = static_cast<uint32_t>(
        std::chrono::system_clock::now().time_since_epoch().count());
    
    // Create plate generator
    AeonTerra::Tools::PlateGenerator plateGen(seed);
    std::cout << "Generating " << plateCount << " tectonic plates..." << std::endl;
    plateGen.generatePlates(plateCount, planetRadius);
    
    // Create tectonic simulation
    AeonTerra::Simulation::PlateTectonics tectonics;
    tectonics.LoadPlatesFromGenerator(plateGen);
    
    // Create heightmap generator
    AeonTerra::Tools::HeightmapGenerator heightmapGen;
    AeonTerra::Tools::HeightmapGenerator::Options options;
    options.resolution = resolution;
    
    // Run simulation and generate frames
    float frameTimeYears = simulationYears / frameCount;
    
    for (int frame = 0; frame <= frameCount; ++frame) {
        std::cout << "Generating frame " << frame << " at " 
                  << std::fixed << std::setprecision(2) 
                  << (frame * frameTimeYears) << " million years..." << std::endl;
        
        // Generate heightmap from current plate state
        std::vector<float> heightmap = heightmapGen.GenerateHeightmap(tectonics, options);
        
        // Save heightmap as PPM image for visualization
        std::string filename = "planet_" + std::to_string(frame) + ".ppm";
        writePPM(filename, heightmap, resolution, resolution);
        
        // Save raw heightmap data
        heightmapGen.SaveHeightmap("heightmap_" + std::to_string(frame) + ".raw", 
                                  heightmap, resolution);
        
        // Advance simulation
        if (frame < frameCount) {
            tectonics.SimulateFrame(frameTimeYears);
        }
    }
    
    std::cout << "Simulation complete! Generated " << (frameCount + 1) << " frames." << std::endl;
    
    return 0;
}