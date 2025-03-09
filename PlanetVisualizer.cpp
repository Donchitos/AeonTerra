#include "Core/Simulation/PlateTectonics.h"
#include "Core/Simulation/GeologicalTimeline.h"
#include "Tools/WorldBuilder/PlateGenerator.h"
#include "Tools/HeightmapGenerator.h"
#include "Tools/TectonicVisualizer.h"
#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <filesystem>

using namespace AeonTerra;

void PrintUsage() {
    std::cout << "AeonTerra Planet Visualizer\n";
    std::cout << "Usage:\n";
    std::cout << "  -p <plateCount>   : Number of tectonic plates (default: 12)\n";
    std::cout << "  -r <radiusKm>     : Planet radius in km (default: 6371)\n";
    std::cout << "  -s <seed>         : Random seed for generation (default: time-based)\n";
    std::cout << "  -t <timeStepMY>   : Time step in millions of years (default: 10)\n";
    std::cout << "  -o <outputFolder> : Output folder for visualizations (default: ./output)\n";
    std::cout << "  -pangea           : Start with a Pangea-like supercontinent\n";
    std::cout << "  -fantasy          : Create a fantasy world instead of Earth-like\n";
    std::cout << "  -h, --help        : Show this help message\n";
}

int main(int argc, char** argv) {
    // Default parameters
    int plateCount = 12;
    float planetRadius = 6371.0f;
    uint32_t seed = static_cast<uint32_t>(std::chrono::system_clock::now().time_since_epoch().count());
    float timeStepMY = 10.0f;  // Million years
    std::string outputFolder = "./output";
    bool startWithPangea = false;
    bool fantasyWorld = false;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-p" && i + 1 < argc) {
            plateCount = std::stoi(argv[++i]);
        }
        else if (arg == "-r" && i + 1 < argc) {
            planetRadius = std::stof(argv[++i]);
        }
        else if (arg == "-s" && i + 1 < argc) {
            seed = static_cast<uint32_t>(std::stoul(argv[++i]));
        }
        else if (arg == "-t" && i + 1 < argc) {
            timeStepMY = std::stof(argv[++i]);
        }
        else if (arg == "-o" && i + 1 < argc) {
            outputFolder = argv[++i];
        }
        else if (arg == "-pangea") {
            startWithPangea = true;
        }
        else if (arg == "-fantasy") {
            fantasyWorld = true;
        }
        else if (arg == "-h" || arg == "--help") {
            PrintUsage();
            return 0;
        }
    }
    
    std::cout << "AeonTerra Planet Visualizer\n";
    std::cout << "----------------------------\n";
    std::cout << "Plate Count: " << plateCount << "\n";
    std::cout << "Planet Radius: " << planetRadius << " km\n";
    std::cout << "Random Seed: " << seed << "\n";
    std::cout << "Time Step: " << timeStepMY << " million years\n";
    std::cout << "Output Folder: " << outputFolder << "\n";
    std::cout << "Start with Pangea: " << (startWithPangea ? "Yes" : "No") << "\n";
    std::cout << "World Type: " << (fantasyWorld ? "Fantasy" : "Earth-like") << "\n";
    std::cout << "----------------------------\n";
    
    // Create output directory
    std::filesystem::path outputPath(outputFolder);
    if (!std::filesystem::exists(outputPath)) {
        std::filesystem::create_directories(outputPath);
    }
    
    // Initialize planet
    std::cout << "Generating plates...\n";
    Tools::PlateGenerator plateGenerator(seed);
    
    if (fantasyWorld) {
        // Set fantasy parameters
        Tools::PlateGenerator::GenerationParams params;
        params.continentalPercent = 0.4f;    // More continental crust
        params.maxFragmentation = 0.8f;      // More fragmented plates
        params.irregularity = 0.7f;          // More irregular boundaries
        plateGenerator.setGenerationParams(params);
    }
    
    // Generate plates
    plateGenerator.generateEarthLikePlates(plateCount, planetRadius, startWithPangea);
    
    // Initialize simulation
    std::cout << "Setting up tectonic simulation...\n";
    Simulation::PlateTectonics tectonics;
    tectonics.LoadPlatesFromGenerator(plateGenerator);
    
    // Initialize timeline
    Simulation::GeologicalTimeline timeline;
    if (fantasyWorld) {
        timeline.InitializeFantasy(seed);
    } else {
        timeline.InitializeEarthLike();
    }
    
    // Initialize visualizer
    Tools::TectonicVisualizer visualizer;
    
    // Save initial state
    std::cout << "Generating initial visualizations...\n";
    
    // Generate heightmap
    Tools::HeightmapGenerator heightmapGen;
    heightmapGen.SetSeed(seed);
    
    Tools::HeightmapGenerator::Options options;
    options.resolution = 2048;
    
    std::cout << "Generating initial heightmap...\n";
    std::vector<float> initialHeightmap = heightmapGen.GenerateHeightmap(tectonics, options);
    
    // Save initial state visualizations
    std::string initialHeightmapFile = outputFolder + "/initial_heightmap.png";
    std::cout << "Saving initial heightmap to " << initialHeightmapFile << "\n";
    heightmapGen.SaveHeightmapPNG(initialHeightmapFile, initialHeightmap, options.resolution);
    
    std::string initialPlateFile = outputFolder + "/initial_plates.png";
    std::cout << "Saving initial plate map to " << initialPlateFile << "\n";
    heightmapGen.SavePlateMap(initialPlateFile, tectonics, options.resolution);
    
    std::string initialReliefFile = outputFolder + "/initial_relief.png";
    std::cout << "Saving initial relief map to " << initialReliefFile << "\n";
    heightmapGen.SaveShadedReliefMap(initialReliefFile, initialHeightmap, options.resolution);
    
    // Generate cross-section
    std::string crossSectionFile = outputFolder + "/cross_section.png";
    std::cout << "Generating cross-section to " << crossSectionFile << "\n";
    visualizer.GenerateCrustMantelSection(tectonics, 0.0f, 0.0f, 3000.0f, crossSectionFile);
    
    // Generate stress map
    std::string stressFile = outputFolder + "/stress_map.png";
    std::cout << "Generating stress map to " << stressFile << "\n";
    visualizer.GenerateStressMap(tectonics, stressFile);
    
    // Simulate tectonic evolution
    // Number of steps to simulate
    const int SIMULATION_STEPS = 10;
    
    std::cout << "\nSimulating tectonic evolution...\n";
    for (int step = 0; step < SIMULATION_STEPS; step++) {
        float progress = static_cast<float>(step) / SIMULATION_STEPS * 100.0f;
        std::cout << "Simulation progress: " << progress << "%" << std::endl;
        
        // Run simulation for this time step
        float yearsToSimulate = timeStepMY * 1000000.0f;
        tectonics.SimulateFrame(yearsToSimulate);
        
        // Save state at this step
        std::string stepFolder = outputFolder + "/step_" + std::to_string(step + 1);
        std::filesystem::create_directories(stepFolder);
        
        // Generate heightmap for current step
        std::vector<float> currentHeightmap = heightmapGen.GenerateHeightmap(tectonics, options);
        
        // Save visualizations
        std::string heightmapFile = stepFolder + "/heightmap.png";
        heightmapGen.SaveHeightmapPNG(heightmapFile, currentHeightmap, options.resolution);
        
        std::string plateFile = stepFolder + "/plates.png";
        heightmapGen.SavePlateMap(plateFile, tectonics, options.resolution);
        
        std::string reliefFile = stepFolder + "/relief.png";
        heightmapGen.SaveShadedReliefMap(reliefFile, currentHeightmap, options.resolution);
        
        std::string crossFile = stepFolder + "/cross_section.png";
        visualizer.GenerateCrustMantelSection(tectonics, 0.0f, 0.0f, 3000.0f, crossFile);
        
        std::string stressFile = stepFolder + "/stress_map.png";
        visualizer.GenerateStressMap(tectonics, stressFile);
    }
    
    std::cout << "\nSimulation complete! Results saved to " << outputFolder << "\n";
    std::cout << "Open the images to see the planet's evolution over time.\n";
    
    return 0;
}