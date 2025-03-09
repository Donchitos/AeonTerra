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
            , oceanDepth(-8000.0f)
            , continentHeight(400.0f)
            , mountainScale(8000.0f)
            , noiseScale(500.0f)
        {}
    };
    
    HeightmapGenerator();
    
    // Set the random seed for noise generation
    void SetSeed(uint32_t seed);
    
    // Generate a heightmap from plate tectonic simulation
    std::vector<float> GenerateHeightmap(
        const AeonTerra::Simulation::PlateTectonics& tectonics,
        const Options& options = Options()
    );
    
    // Save raw heightmap to binary file
    bool SaveHeightmap(
        const std::string& filename, 
        const std::vector<float>& heightmap, 
        int resolution
    );
    
    // Save heightmap as a colored PNG image
    bool SaveHeightmapPNG(
        const std::string& filename,
        const std::vector<float>& heightmap,
        int resolution
    );
    
    // Save a map showing plate boundaries
    bool SavePlateMap(
        const std::string& filename,
        const AeonTerra::Simulation::PlateTectonics& tectonics,
        int resolution
    );
    
    // Save a before/after comparison image
    bool SaveComparisonImage(
        const std::string& filename,
        const std::vector<float>& beforeHeightmap,
        const std::vector<float>& afterHeightmap,
        int resolution
    );
    
    // Save a 3D shaded relief map with lighting
    bool SaveShadedReliefMap(
        const std::string& filename,
        const std::vector<float>& heightmap,
        int resolution
    );
    
private:
    // Calculate base height at a specific point on the sphere based on plate properties
    float CalculateBaseHeight(
        const AeonTerra::Math::Vector3& point, 
        const AeonTerra::Simulation::PlateTectonics& tectonics,
        const Options& options
    );
    
    // Calculate height contribution from plate boundaries
    float CalculateBoundaryHeight(
        const AeonTerra::Math::Vector3& point,
        const std::vector<AeonTerra::Simulation::PlateTectonics::BoundaryPoint>& boundaries,
        const Options& options
    );
    
    // Generate fractal noise for terrain variation
    float CalculateFractalNoise(
        const AeonTerra::Math::Vector3& point, 
        const std::vector<float>& noiseGrid, 
        int gridSize,
        float scale
    );
    
    // Normalize heightmap to specified range
    void NormalizeHeightmap(
        std::vector<float>& heightmap,
        float minValue,
        float maxValue
    );
    
    // Random seed for noise generation
    uint32_t m_seed;
};

} // namespace Tools
} // namespace AeonTerra