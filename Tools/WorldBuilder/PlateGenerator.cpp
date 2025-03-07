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

// Rest of existing implementation...