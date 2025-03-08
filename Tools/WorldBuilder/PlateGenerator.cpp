#include "Tools/WorldBuilder/PlateGenerator.h"
#include "Core/CoreExports.h"
#include "Core/Mathematics/GeoMath.h"
#include "Core/Mathematics/Vector3.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <numbers>
#include <random>
#include <stdexcept>
#include <vector>

namespace AeonTerra::Tools {

AeonTerra::Tools::PlateGenerator::PlateGenerator(uint32_t seed)
    : m_engine(seed),
      m_plateCount(0),
      m_planetRadius(6371.0f),
      m_plates() {}

void PlateGenerator::generatePlates(int plateCount, float planetRadius) {
    if(plateCount < 1 || planetRadius <= 0) {
        throw std::invalid_argument("Invalid plate generation parameters");
    }

    std::vector<Vector3> points;
    
    // Generate random points on sphere surface
    for(int i = 0; i < plateCount; ++i) {
        points.push_back(GeoMath::randomUnitVector(m_engine));
    }

    // Create Voronoi cells
    m_plates.clear();
    for(const auto& point : points) {
        Plate plate;
        plate.center = point * planetRadius;
        plate.velocity = GeoMath::randomUnitVector(m_engine) * 0.5f;
        plate.buoyancy = std::uniform_real_distribution<float>(0.7f, 1.3f)(m_engine);
        m_plates.push_back(plate);
    }

    CreateVoronoiCells(m_plates, planetRadius);
}

void PlateGenerator::CreateVoronoiCells(std::vector<Plate>& plates, float sphereRadius) {
    // Generate Voronoi seed points
    std::vector<SphericalCoord> seeds;
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
            SphericalCoord testPoint{
                sphereRadius,
                thetaDist(m_engine),
                phiDist(m_engine)
            };
            
            // Find nearest seed point
            auto nearest = std::min_element(seeds.begin(), seeds.end(),
                [&](const SphericalCoord& a, const SphericalCoord& b) {
                    return GeoMath::GreatCircleDistance(testPoint, a)
                         < GeoMath::GreatCircleDistance(testPoint, b);
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

} // namespace AeonTerra::Tools
#include "PlateGenerator.h"
#include "../../Core/Mathematics/GeoMath.h"
#include <random>
#include <cmath>
#include <algorithm>

namespace AeonTerra {
namespace Tools {

using namespace Math;
using namespace Simulation;

std::vector<PlateTectonics::Plate> PlateGenerator::GeneratePlates(int count, float sphereRadius) {
    std::vector<PlateTectonics::Plate> plates(count);
    
    // Generate random seed points on sphere using spherical coordinates
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> thetaDist(0.0f, M_PI);
    std::uniform_real_distribution<float> phiDist(0.0f, 2*M_PI);
    
    CreateVoronoiCells(plates, sphereRadius);
    AssignPlateProperties(plates);
    
    return plates;
}

void PlateGenerator::CreateVoronoiCells(std::vector<PlateTectonics::Plate>& plates, float sphereRadius) {
    // Generate Voronoi seed points
    std::vector<SphericalCoord> seeds;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> thetaDist(0.0f, M_PI);
    std::uniform_real_distribution<float> phiDist(0.0f, 2*M_PI);
    
    for(int i=0; i<plates.size(); i++) {
        seeds.push_back({sphereRadius, thetaDist(gen), phiDist(gen)});
        plates[i].id = i;
    }
    
    // Simplified Voronoi generation using nearest neighbor search
    const int boundaryPoints = 100;
    for(auto& plate : plates) {
        plate.boundaries.clear();
        for(int i=0; i<boundaryPoints; i++) {
            SphericalCoord testPoint{
                sphereRadius,
                thetaDist(gen),
                phiDist(gen)
            };
            
            // Find nearest seed point
            auto nearest = std::min_element(seeds.begin(), seeds.end(),
                [&](const SphericalCoord& a, const SphericalCoord& b) {
                    return GeoMath::GreatCircleDistance(testPoint, a) 
                         < GeoMath::GreatCircleDistance(testPoint, b);
                });
            
            if(plate.id == std::distance(seeds.begin(), nearest)) {
                plate.boundaries.push_back(testPoint);
            }
        }
    }
}

void PlateGenerator::AssignPlateProperties(std::vector<PlateTectonics::Plate>& plates) {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Distribution for realistic plate properties
    std::uniform_real_distribution<float> angVelDist(-0.01f, 0.01f);  // radians/year
    std::uniform_real_distribution<float> thicknessDist(30.0f, 100.0f); // km
    std::uniform_real_distribution<float> densityDist(2700.0f, 3300.0f); // kg/m³
    
    for(auto& plate : plates) {
        plate.angularVelocity = angVelDist(gen);
        plate.thickness = thicknessDist(gen);
        plate.density = densityDist(gen);
    }
}
// Rest of existing implementation...