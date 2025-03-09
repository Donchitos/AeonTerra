#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <filesystem>
#include "Core/Simulation/PlateTectonics.h"
#include "Tools/WorldBuilder/PlateGenerator.h"
#include "Tools/HeightmapGenerator.h"

// Visualization configuration
struct VisualizationConfig {
    int resolution = 1024;          // Image resolution
    int plateCount = 15;            // Number of tectonic plates
    float planetRadius = 6371.0f;   // Planet radius in km (Earth-like)
    float timeStep = 1000000.0f;    // Simulation time step in years
    int frames = 20;                // Number of frames to simulate
    std::string outputDir = "./output"; // Output directory for images
    bool showPlates = true;         // Generate plate boundary maps
    bool showHeightmap = true;      // Generate height maps
    bool showRelief = true;         // Generate relief maps
    bool createAnimation = true;    // Create animation frames
    bool useRealisticPlates = true; // Use the new realistic plate generation
    bool startWithPangea = false;   // Start with a Pangea-like supercontinent
};

void printSimulationInfo(const AeonTerra::Simulation::PlateTectonics& tectonics) {
    const auto& plates = tectonics.GetPlates();
    const auto& boundaries = tectonics.GetConvergentBoundaries();
    
    std::cout << "=== Simulation Information ===" << std::endl;
    std::cout << "Time: " << tectonics.GetSimulationTime() / 1000000.0 << " million years" << std::endl;
    std::cout << "Plates: " << plates.size() << std::endl;
    std::cout << "Convergent Boundaries: " << boundaries.size() << std::endl;
    
    // Show plate details
    std::cout << "\nPlate Details:" << std::endl;
    std::cout << "ID\tDensity\tThickness\tAngular Velocity" << std::endl;
    for (const auto& plate : plates) {
        std::cout << plate.id << "\t" 
                  << plate.density << "\t" 
                  << plate.thickness << "\t\t" 
                  << plate.angularVelocity << std::endl;
    }
    
    // Show high mountains and deep trenches
    if (!boundaries.empty()) {
        float maxHeight = boundaries[0].height;
        float minHeight = boundaries[0].height;
        int maxIdx = 0, minIdx = 0;
        
        for (size_t i = 1; i < boundaries.size(); i++) {
            if (boundaries[i].height > maxHeight) {
                maxHeight = boundaries[i].height;
                maxIdx = i;
            }
            if (boundaries[i].height < minHeight) {
                minHeight = boundaries[i].height;
                minIdx = i;
            }
        }
        
        std::cout << "\nHighest mountain: " << maxHeight << " meters" << std::endl;
        std::cout << "Deepest trench: " << minHeight << " meters" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    VisualizationConfig config;
    
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(config.outputDir);
    
    std::cout << "AeonTerra Planet Visualizer" << std::endl;
    std::cout << "==========================" << std::endl;
    
    // Initialize plate generator with random seed
    std::random_device rd;
    AeonTerra::Tools::PlateGenerator plateGen(rd());
    
    // Set custom parameters for realistic plate generation
    AeonTerra::Tools::PlateGenerator::GenerationParams genParams;
    genParams.continentalPercent = 0.35f;  // 35% continental crust (like Earth)
    genParams.maxFragmentation = 0.7f;     // Degree of fragmentation
    genParams.irregularity = 0.65f;        // Irregularity of plate shapes
    plateGen.setGenerationParams(genParams);
    
    // Generate initial plates using the improved realistic method
    std::cout << "Generating " << config.plateCount << " tectonic plates..." << std::endl;
    if (config.useRealisticPlates) {
        plateGen.generateEarthLikePlates(config.plateCount, config.planetRadius, config.startWithPangea);
        std::cout << "Using realistic Earth-like plate generation" << std::endl;
    } else {
        plateGen.generatePlates(config.plateCount, config.planetRadius);
        std::cout << "Using basic plate generation" << std::endl;
    }
    
    // Create tectonic simulation
    AeonTerra::Simulation::PlateTectonics tectonics;
    tectonics.LoadPlatesFromGenerator(plateGen);
    
    // Configure tectonic parameters
    tectonics.SetMantleViscosity(1e21f);  // Earth-like mantle viscosity
    tectonics.SetSeafloorSpreadingRate(50.0f);    // 50 mm/year (typical Earth rate)
    tectonics.SetSubductionRate(80.0f);           // 80 mm/year
    tectonics.SetContinentalCollisionRate(30.0f); // 30 mm/year
    
    // Initialize heightmap generator
    AeonTerra::Tools::HeightmapGenerator heightmapGen;
    heightmapGen.SetSeed(rd());
    
    // Configure heightmap options
    AeonTerra::Tools::HeightmapGenerator::Options heightmapOptions;
    heightmapOptions.resolution = config.resolution;
    heightmapOptions.oceanDepth = -8000.0f;
    heightmapOptions.continentHeight = 500.0f;
    heightmapOptions.mountainScale = 8000.0f;
    
    // Generate initial heightmap
    std::vector<float> initialHeightmap = heightmapGen.GenerateHeightmap(tectonics, heightmapOptions);
    
    // Save initial state
    std::string initialHeightmapFile = config.outputDir + "/initial_heightmap.ppm";
    std::string initialPlatesFile = config.outputDir + "/initial_plates.ppm";
    std::string initialReliefFile = config.outputDir + "/initial_relief.ppm";
    
    if (config.showHeightmap) {
        std::cout << "Saving initial heightmap..." << std::endl;
        heightmapGen.SaveHeightmapPNG(initialHeightmapFile, initialHeightmap, config.resolution);
    }
    
    if (config.showPlates) {
        std::cout << "Saving initial plate map..." << std::endl;
        heightmapGen.SavePlateMap(initialPlatesFile, tectonics, config.resolution);
    }
    
    if (config.showRelief) {
        std::cout << "Saving initial relief map..." << std::endl;
        heightmapGen.SaveShadedReliefMap(initialReliefFile, initialHeightmap, config.resolution);
    }
    
    // Run simulation over multiple time steps
    if (config.frames > 1) {
        std::cout << "\nRunning plate tectonic simulation..." << std::endl;
        std::cout << "Simulating " << config.frames << " frames with " 
                  << config.timeStep / 1000000.0 << " million years per step" << std::endl;
        
        // Store a copy of the initial heightmap for final comparison
        std::vector<float> veryInitialHeightmap = initialHeightmap;
        
        for (int frame = 1; frame <= config.frames; frame++) {
            std::cout << "\nSimulating frame " << frame << "..." << std::endl;
            
            // Run simulation step
            tectonics.SimulateFrame(config.timeStep);
            
            // Print current simulation state
            printSimulationInfo(tectonics);
            
            // Generate new heightmap
            std::vector<float> currentHeightmap = heightmapGen.GenerateHeightmap(tectonics, heightmapOptions);
            
            if (config.createAnimation) {
                // Save heightmap for this frame
                std::string frameHeightmapFile = config.outputDir + "/frame_" + 
                                              std::to_string(frame) + "_heightmap.ppm";
                std::string framePlatesFile = config.outputDir + "/frame_" + 
                                           std::to_string(frame) + "_plates.ppm";
                std::string frameReliefFile = config.outputDir + "/frame_" + 
                                           std::to_string(frame) + "_relief.ppm";
                
                if (config.showHeightmap) {
                    heightmapGen.SaveHeightmapPNG(frameHeightmapFile, currentHeightmap, config.resolution);
                }
                
                if (config.showPlates) {
                    heightmapGen.SavePlateMap(framePlatesFile, tectonics, config.resolution);
                }
                
                if (config.showRelief) {
                    heightmapGen.SaveShadedReliefMap(frameReliefFile, currentHeightmap, config.resolution);
                }
                
                // Save comparison with previous frame
                if (frame > 1) {
                    std::string comparisonFile = config.outputDir + "/comparison_" + 
                                              std::to_string(frame-1) + "_to_" + 
                                              std::to_string(frame) + ".ppm";
                    heightmapGen.SaveComparisonImage(comparisonFile, initialHeightmap, 
                                                   currentHeightmap, config.resolution);
                }
            }
            
            // Use current heightmap as reference for next comparison
            initialHeightmap = currentHeightmap;
            
            std::cout << "Frame " << frame << " complete." << std::endl;
        }
        
        // Save final comparison between very first frame and last frame
        std::string finalComparisonFile = config.outputDir + "/final_comparison.ppm";
        heightmapGen.SaveComparisonImage(finalComparisonFile, veryInitialHeightmap, initialHeightmap, config.resolution);
    }
    
    std::cout << "\nSimulation complete. Output files saved to: " << config.outputDir << std::endl;
    std::cout << "To view the results, open the PPM files with an image viewer that supports this format." << std::endl;
    std::cout << "Alternatively, convert them to PNG using ImageMagick or similar tools." << std::endl;
    
    return 0;
}