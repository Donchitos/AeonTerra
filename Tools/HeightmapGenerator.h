#pragma once
#include "../Core/CoreExports.h"
#include "../Core/Simulation/PlateTectonics.h"
#include "../Core/Mathematics/Vector3.h"
#include <vector>
#include <memory>
#include <string>
#include <limits>

namespace AeonTerra {
namespace Tools {

class AEONTERRA_API HeightmapGenerator {
public:
    // Configuration options for heightmap generation
    struct Options {
        int resolution;           // Resolution of heightmap (width/height)
        float oceanDepth;         // Base depth of oceans (negative height)
        float continentHeight;    // Base height of continents
        float mountainScale;      // Scale factor for mountain heights
        float noiseScale;         // Scale of noise to add to heightmap
        
        // Default constructor with reasonable values
        Options() 
            : resolution(1024)
            , oceanDepth(-4000.0f)
            , continentHeight(400.0f)
            , mountainScale(8000.0f)
            , noiseScale(500.0f)
        {}
    };
    
    HeightmapGenerator();
    
    // Generate a heightmap from plate tectonic simulation
    std::vector<float> GenerateHeightmap(
        const AeonTerra::Simulation::PlateTectonics& tectonics,
        const Options& options = Options()
    );
    
    // Save heightmap to file
    bool SaveHeightmap(
        const std::string& filename, 
        const std::vector<float>& heightmap, 
        int resolution
    );
    
private:
    // Calculate height at a specific point on the sphere
    float CalculateHeight(
        const AeonTerra::Math::Vector3& point, 
        const AeonTerra::Simulation::PlateTectonics& tectonics,
        const Options& options
    );
    
    // Calculate distance to nearest plate boundary
    float DistanceToPlateBoundary(
        const AeonTerra::Math::Vector3& point,
        const AeonTerra::Simulation::PlateTectonics& tectonics
    );
    
    // Generate simple noise for terrain variation
    float GenerateNoise(float x, float y, float scale);
};

} // namespace Tools
} // namespace AeonTerra