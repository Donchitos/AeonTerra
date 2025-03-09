#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>
#include <thread>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>

// Define M_PI if not available (for MSVC compatibility)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Core/Simulation/PlateTectonics.h"
#include "Tools/WorldBuilder/PlateGenerator.h"
#include "Tools/HeightmapGenerator.h"

namespace fs = std::filesystem;

// Helper function to create a timestamped filename
std::string CreateTimestampedFilename(const std::string& prefix, const std::string& extension, int timestep) {
    std::ostringstream oss;
    oss << prefix << "_" << std::setw(6) << std::setfill('0') << timestep << "." << extension;
    return oss.str();
}

// Attempt to create a video from the generated images
bool CreateAnimation(const std::string& outputDir, const std::string& prefix, int startYear, int endYear, int interval) {
    // Check if ffmpeg is available
    std::string checkCommand;
    #ifdef _WIN32
    checkCommand = "where ffmpeg > nul 2>&1";
    #else
    checkCommand = "which ffmpeg > /dev/null 2>&1";
    #endif
    
    int checkResult = system(checkCommand.c_str());
    if (checkResult != 0) {
        std::cout << "FFmpeg not found. Cannot create animation.\n";
        std::cout << "To create an animation, install FFmpeg and run:\n";
        std::cout << "ffmpeg -framerate 5 -i " << outputDir << "/" << prefix << "_%06d.ppm -c:v libx264 -pix_fmt yuv420p " << outputDir << "/animation.mp4\n";
        return false;
    }
    
    // Construct FFmpeg command
    std::string cmd;
    #ifdef _WIN32
    cmd = "ffmpeg -y -framerate 5 -start_number " + std::to_string(startYear) + " -i \"" + outputDir + "/" + prefix + "_%06d.ppm\" -c:v libx264 -pix_fmt yuv420p \"" + outputDir + "/animation_" + prefix + ".mp4\"";
    #else
    cmd = "ffmpeg -y -framerate 5 -start_number " + std::to_string(startYear) + " -i '" + outputDir + "/" + prefix + "_%06d.ppm' -c:v libx264 -pix_fmt yuv420p '" + outputDir + "/animation_" + prefix + ".mp4'";
    #endif
    
    std::cout << "Creating animation with command: " << cmd << std::endl;
    int result = system(cmd.c_str());
    
    if (result != 0) {
        std::cout << "Error creating animation. FFmpeg returned error code: " << result << std::endl;
        return false;
    }
    
    std::cout << "Animation created: " << outputDir << "/animation_" << prefix << ".mp4" << std::endl;
    return true;
}

int main(int argc, char* argv[]) {
    std::cout << "AeonTerra Planet Simulator\n";
    std::cout << "==========================\n\n";
    
    // Create output directory
    std::string outputDir = "simulation_output";
    if (!fs::exists(outputDir)) {
        fs::create_directory(outputDir);
    }
    
    // Parse command line args or use defaults
    int plateCount = 12;               // Default number of plates
    float planetRadius = 6371.0f;      // Earth radius in km
    int simulationYears = 50;          // Years to simulate
    int snapshotInterval = 5;          // Take snapshot every N years
    int resolution = 800;              // Resolution of output images
    uint32_t seed = 12345;             // Random seed
    bool createAnimation = true;       // Whether to create animation
    
    // Parse command line arguments if provided
    if (argc > 1) plateCount = std::stoi(argv[1]);
    if (argc > 2) simulationYears = std::stoi(argv[2]);
    if (argc > 3) snapshotInterval = std::stoi(argv[3]);
    if (argc > 4) resolution = std::stoi(argv[4]);
    if (argc > 5) seed = std::stoul(argv[5]);
    if (argc > 6) createAnimation = (std::stoi(argv[6]) != 0);
    
    std::cout << "Configuration:\n";
    std::cout << "  Plate Count: " << plateCount << "\n";
    std::cout << "  Planet Radius: " << planetRadius << " km\n";
    std::cout << "  Simulation Years: " << simulationYears << "\n";
    std::cout << "  Snapshot Interval: " << snapshotInterval << " years\n";
    std::cout << "  Resolution: " << resolution << "x" << resolution << "\n";
    std::cout << "  Random Seed: " << seed << "\n";
    std::cout << "  Create Animation: " << (createAnimation ? "Yes" : "No") << "\n\n";

    // Generate initial plates
    std::cout << "Generating initial plates...\n";
    AeonTerra::Tools::PlateGenerator plateGenerator(seed);
    plateGenerator.generatePlates(plateCount, planetRadius);
    
    // Initialize plate tectonics simulation
    std::cout << "Initializing plate tectonics simulation...\n";
    AeonTerra::Simulation::PlateTectonics tectonics;
    tectonics.LoadPlatesFromGenerator(plateGenerator);
    tectonics.SetMantleViscosity(1e21f); // Set mantle viscosity - affects plate speed
    
    // Initialize heightmap generator with the same seed
    AeonTerra::Tools::HeightmapGenerator heightmapGen;
    heightmapGen.SetSeed(seed);
    
    // Configure heightmap options
    AeonTerra::Tools::HeightmapGenerator::Options heightmapOptions;
    heightmapOptions.resolution = resolution;
    heightmapOptions.oceanDepth = -8000.0f;
    heightmapOptions.continentHeight = 400.0f;
    heightmapOptions.mountainScale = 8000.0f;
    heightmapOptions.noiseScale = 500.0f;
    
    // Generate and save initial state
    std::cout << "Generating initial heightmap...\n";
    auto initialHeightmap = heightmapGen.GenerateHeightmap(tectonics, heightmapOptions);
    
    std::string heightmapFilename = outputDir + "/" + CreateTimestampedFilename("heightmap", "ppm", 0);
    std::string shadedReliefFilename = outputDir + "/" + CreateTimestampedFilename("relief", "ppm", 0);
    std::string plateMapFilename = outputDir + "/" + CreateTimestampedFilename("platemap", "ppm", 0);
    
    std::cout << "Saving initial visualization to " << heightmapFilename << "...\n";
    heightmapGen.SaveHeightmapPNG(heightmapFilename, initialHeightmap, resolution);
    heightmapGen.SaveShadedReliefMap(shadedReliefFilename, initialHeightmap, resolution);
    heightmapGen.SavePlateMap(plateMapFilename, tectonics, resolution);
    
    // Store all heightmaps for later comparison and animation
    std::vector<std::vector<float>> heightmapHistory;
    heightmapHistory.push_back(initialHeightmap);
    
    // Main simulation loop
    std::cout << "\nStarting simulation...\n";
    for (int year = 1; year <= simulationYears; ++year) {
        std::cout << "Simulating year " << year << "/" << simulationYears << "... ";
        
        // Advance the simulation by 1 year
        tectonics.SimulateFrame(1.0f);
        
        // Take snapshots at specified intervals
        if (year % snapshotInterval == 0) {
            std::cout << "Taking snapshot... ";
            
            // Generate new heightmap based on current state
            auto currentHeightmap = heightmapGen.GenerateHeightmap(tectonics, heightmapOptions);
            heightmapHistory.push_back(currentHeightmap);
            
            // Save visualization files
            heightmapFilename = outputDir + "/" + CreateTimestampedFilename("heightmap", "ppm", year);
            shadedReliefFilename = outputDir + "/" + CreateTimestampedFilename("relief", "ppm", year);
            plateMapFilename = outputDir + "/" + CreateTimestampedFilename("platemap", "ppm", year);
            
            heightmapGen.SaveHeightmapPNG(heightmapFilename, currentHeightmap, resolution);
            heightmapGen.SaveShadedReliefMap(shadedReliefFilename, currentHeightmap, resolution);
            heightmapGen.SavePlateMap(plateMapFilename, tectonics, resolution);
            
            // Also save a before/after comparison image
            std::string comparisonFilename = outputDir + "/" + CreateTimestampedFilename("comparison", "ppm", year);
            heightmapGen.SaveComparisonImage(comparisonFilename, initialHeightmap, currentHeightmap, resolution);
            
            std::cout << "Done.\n";
        } else {
            std::cout << "Done.\n";
        }
    }
    
    // Create animations if requested
    if (createAnimation) {
        std::cout << "\nCreating animations...\n";
        CreateAnimation(outputDir, "heightmap", 0, simulationYears, snapshotInterval);
        CreateAnimation(outputDir, "relief", 0, simulationYears, snapshotInterval);
        CreateAnimation(outputDir, "platemap", 0, simulationYears, snapshotInterval);
    }
    
    std::cout << "\nSimulation complete!\n";
    std::cout << "Output files saved to: " << fs::absolute(outputDir).string() << "\n";
    
    // Display instructions for viewing
    std::cout << "\nTo view the results:\n";
    std::cout << "1. Open the PPM files with an image viewer (GIMP, IrfanView, etc.)\n";
    std::cout << "2. Look at the comparison images to see before/after changes\n";
    if (createAnimation) {
        std::cout << "3. Watch the animations to see plate movement over time\n";
    } else {
        std::cout << "3. To create animations manually, install FFmpeg and run:\n";
        std::cout << "   ffmpeg -framerate 5 -i " << outputDir << "/heightmap_%06d.ppm -c:v libx264 -pix_fmt yuv420p " << outputDir << "/animation.mp4\n";
    }
    
    return 0;
}