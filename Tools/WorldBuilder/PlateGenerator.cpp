#include "PlateGenerator.h"
#include "Core/Mathematics/GeoMath.h"
#include "Core/Mathematics/Vector3.h"
#include <random>
#include <cmath>
#include <algorithm>

namespace AeonTerra {
namespace Tools {

PlateGenerator::PlateGenerator(uint32_t seed)
    : m_engine(seed),
      m_plateCount(0),
      m_planetRadius(6371.0f),
      m_plates() {}

void PlateGenerator::generatePlates(int plateCount, float planetRadius) {
    if(plateCount < 1 || planetRadius <= 0) {
        throw std::invalid_argument("Invalid plate generation parameters");
    }

    m_plateCount = plateCount;
    m_planetRadius = planetRadius;
    
    std::vector<Math::Vector3> points;
    
    // Generate random points on sphere surface
    for(int i = 0; i < plateCount; ++i) {
        points.push_back(Math::GeoMath::randomUnitVector(m_engine));
    }

    // Create Voronoi cells
    m_plates.clear();
    for(int i = 0; i < points.size(); ++i) {
        Plate plate;
        plate.id = i;
        plate.center = points[i] * planetRadius;
        plate.velocity = Math::GeoMath::randomUnitVector(m_engine) * 0.5f;
        plate.buoyancy = std::uniform_real_distribution<float>(0.7f, 1.3f)(m_engine);
        
        // Add plate thickness and density for compatibility with PlateTectonics
        std::uniform_real_distribution<float> thicknessDist(30.0f, 100.0f); // km
        std::uniform_real_distribution<float> densityDist(2700.0f, 3300.0f); // kg/m³
        plate.thickness = thicknessDist(m_engine);
        plate.density = densityDist(m_engine);
        
        m_plates.push_back(plate);
    }

    CreateVoronoiCells(m_plates, planetRadius);
}

void PlateGenerator::CreateVoronoiCells(std::vector<Plate>& plates, float sphereRadius) {
    // Generate Voronoi seed points
    std::vector<Math::SphericalCoord> seeds;
    std::uniform_real_distribution<float> thetaDist(0.0f, static_cast<float>(M_PI));
    std::uniform_real_distribution<float> phiDist(0.0f, static_cast<float>(2*M_PI));
    
    for(auto& plate : plates) {
        seeds.push_back({
            sphereRadius,
            thetaDist(m_engine),
            phiDist(m_engine)
        });
    }
    
    // Simplified Voronoi generation using nearest neighbor search
    const int boundaryPoints = 100;
    for(auto& plate : plates) {
        plate.boundaries.clear();
        for(int i=0; i<boundaryPoints; i++) {
            Math::SphericalCoord testPoint{
                sphereRadius,
                thetaDist(m_engine),
                phiDist(m_engine)
            };
            
            // Find nearest seed point
            auto nearest = std::min_element(seeds.begin(), seeds.end(),
                [&](const Math::SphericalCoord& a, const Math::SphericalCoord& b) {
                    return Math::GeoMath::GreatCircleDistance(testPoint, a)
                         < Math::GeoMath::GreatCircleDistance(testPoint, b);
                });
            
            if(plate.id == static_cast<int>(std::distance(seeds.begin(), nearest))) {
                plate.boundaries.push_back(testPoint);
            }
        }
    }
}

const std::vector<Plate>& PlateGenerator::getPlates() const {
    return m_plates;
}

} // namespace Tools
} // namespace AeonTerra