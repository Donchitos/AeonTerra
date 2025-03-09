#pragma once

#include "../../Core/CoreExports.h"
#include "../../Core/Mathematics/Vector3.h"
#include "../../Core/Mathematics/GeoMath.h"
#include "../../Core/Simulation/PlateTectonics.h"
#include <random>
#include <vector>
#include <stdexcept>
#include <array>
#include <unordered_map>

namespace AeonTerra {
namespace Tools {

// Enhanced plate structure with more realistic properties
struct Plate {
    int id;
    Math::Vector3 center;
    Math::Vector3 velocity;
    float buoyancy;
    float thickness;
    float density;
    bool isOceanic;
    float age;        // Age in million years
    std::vector<Math::SphericalCoord> boundaries;
    std::vector<int> neighborIds;  // IDs of adjacent plates
};

// Fault types for more realistic plate boundaries
enum class FaultType {
    Ridge,      // Divergent - new crust forms
    Subduction, // Convergent - crust destroyed
    Transform,  // Sliding boundary
    Collision,  // Continental collision
    Rift        // Continental splitting
};

// Fracture system for realistic plate generation
struct Fracture {
    Math::Vector3 start;
    Math::Vector3 end;
    float width;
    float depth;
};

class PlateGenerator {
public:
    explicit PlateGenerator(uint32_t seed = std::random_device{}());
    
    // Basic plate generation (for compatibility)
    void generatePlates(int plateCount, float planetRadius);
    
    // New realistic plate generation
    void generateEarthLikePlates(int plateCount, float planetRadius, bool startWithPangea = false);
    
    // Getter for plates
    const std::vector<Plate>& getPlates() const;

    // Convert to PlateTectonics::Plate format
    std::vector<Simulation::PlateTectonics::Plate> convertToTectonicPlates() const;
    
    // Parameters for plate generation
    struct GenerationParams {
        float continentalPercent = 0.3f;    // Percentage of continental crust
        float maxFragmentation = 0.7f;      // Degree of plate fragmentation (0-1)
        float irregularity = 0.65f;         // How irregular plate shapes are (0-1)
        float hotspotDensity = 0.0002f;     // Density of mantle hotspots
        bool includeHotspots = true;        // Whether to include hotspots
    };
    
    // Set generation parameters
    void setGenerationParams(const GenerationParams& params) {
        m_params = params;
    }
    
    // Get current parameters
    const GenerationParams& getGenerationParams() const {
        return m_params;
    }

private:
    // New methods for realistic plate generation
    void createFractureNetwork(float planetRadius);
    void growPlatesFromFractures(int plateCount);
    void assignPlateProperties(bool pangea);
    void createWeakZones();
    void balanceCrustTypes(float continentalPercent);
    void createHotspots();
    
    // Helper methods
    Math::Vector3 calculateInitialVelocity(const Math::Vector3& position, bool isContinental);
    void CreateVoronoiCells(std::vector<Plate>& plates, float sphereRadius);
    float getRandomContinentalThickness();
    float getRandomOceanicThickness();
    float getRandomContinentalDensity();
    float getRandomOceanicDensity();
    
    
    // Improved methods for plate formation
    void initializeCrackedCrust(float planetRadius);
    void applyFractalPerturbations();
    
    // Engine and data
    std::mt19937 m_engine;
    int m_plateCount;
    float m_planetRadius;
    std::vector<Plate> m_plates;
    std::vector<Fracture> m_fractures;
    std::vector<Math::Vector3> m_hotspots;
    GenerationParams m_params;

    private:
    // Add these two helper method declarations
    float distancePointToLineSegment(
        const Math::Vector3& point, 
        const Math::Vector3& lineStart, 
        const Math::Vector3& lineEnd);
        
    bool lineSegmentsIntersect(
        const Math::Vector3& p1, const Math::Vector3& p2,
        const Math::Vector3& p3, const Math::Vector3& p4);
    
};

} // namespace Tools
} // namespace AeonTerra