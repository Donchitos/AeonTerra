#pragma once

#include "Core/CoreExports.h"
#include "Core/Simulation/PlateTectonics.h"
#include "Tools/HeightmapGenerator.h"
#include <string>

namespace AeonTerra {
namespace Tools {

class AEONTERRA_API TectonicVisualizer {
public:
    // Creates time-lapse visualization of plate movement
    void GeneratePlateTectonicAnimation(
        const Simulation::PlateTectonics& tectonics,
        float timeStart, // in million years ago
        float timeEnd,   // in million years ago
        float timeStep,  // in million years
        const std::string& outputDir
    );
    
    // Show cross-section of planetary interior
    void GenerateCrustMantelSection(
        const Simulation::PlateTectonics& tectonics,
        float latitude,  // Cross-section latitude
        float longitude, // Cross-section longitude
        float depth,     // in km
        const std::string& outputFile
    );
    
    // Creates heat map of tectonic stress
    void GenerateStressMap(
        const Simulation::PlateTectonics& tectonics,
        const std::string& outputFile
    );
};

} // namespace Tools
} // namespace AeonTerra