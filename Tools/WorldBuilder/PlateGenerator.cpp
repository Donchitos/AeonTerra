#include "PlateGenerator.h"
#include "../../Core/Mathematics/GeoMath.h"
#include "../../Core/Mathematics/Vector3.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <queue>
#include <unordered_set>

// Define M_PI if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace AeonTerra {
namespace Tools {

PlateGenerator::PlateGenerator(uint32_t seed)
    : m_engine(seed),
      m_plateCount(0),
      m_planetRadius(6371.0f),
      m_plates() {
    // Initialize default generation parameters
    m_params.continentalPercent = 0.3f;
    m_params.maxFragmentation = 0.7f;
    m_params.irregularity = 0.65f;
    m_params.hotspotDensity = 0.0002f;
    m_params.includeHotspots = true;
}

// Original plate generation (kept for compatibility)
void PlateGenerator::generatePlates(int plateCount, float planetRadius) {
    m_planetRadius = planetRadius; 
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
        plate.isOceanic = (plate.density > 3000.0f);
        plate.age = std::uniform_real_distribution<float>(0.0f, 200.0f)(m_engine);
        
        m_plates.push_back(plate);
    }

    CreateVoronoiCells(m_plates, planetRadius);
}

// New Earth-like plate generation
void PlateGenerator::generateEarthLikePlates(int plateCount, float planetRadius, bool startWithPangea) {
    m_planetRadius = planetRadius;
    m_plateCount = plateCount;
    m_plates.clear();
    m_fractures.clear();
    m_hotspots.clear();
    
    // Step 1: Create fracture network that will form plate boundaries
    createFractureNetwork(planetRadius);
    
    // Step 2: Grow plates from fractures to form realistic plate shapes
    growPlatesFromFractures(plateCount);
    
    // Step 3: Assign realistic properties to plates
    assignPlateProperties(startWithPangea);
    
    // Step 4: Create weak zones (potential future boundaries)
    createWeakZones();
    
    // Step 5: Balance continental and oceanic crust distribution
    balanceCrustTypes(m_params.continentalPercent);
    
    // Step 6: Add hotspots if enabled
    if (m_params.includeHotspots) {
        createHotspots();
    }
}

void PlateGenerator::createFractureNetwork(float planetRadius) {
    // Start with a small number of initial cracks
    int initialCracks = static_cast<int>(m_plateCount * 0.5f);
    
    // Create initial "seed" fractures
    for (int i = 0; i < initialCracks; i++) {
        Math::Vector3 startPoint = Math::GeoMath::randomUnitVector(m_engine);
        Math::Vector3 randomDir = Math::GeoMath::randomUnitVector(m_engine);
        Math::Vector3 perpDir = startPoint.cross(randomDir).normalized();
        
        // Initial fracture length 
        float fracLength = std::uniform_real_distribution<float>(0.2f, 0.5f)(m_engine);
        Math::Vector3 endPoint = (startPoint * std::cos(fracLength) + 
                                perpDir * std::sin(fracLength)).normalized();
        
        // Create initial fracture
        m_fractures.push_back({
            startPoint * planetRadius,
            endPoint * planetRadius, 
            std::uniform_real_distribution<float>(0.02f, 0.08f)(m_engine),  // Width
            std::uniform_real_distribution<float>(15.0f, 40.0f)(m_engine)   // Depth
        });
    }
    
    // Propagate fractures using a stress-based system
    const int PROPAGATION_STEPS = m_plateCount * 2;
    
    for (int step = 0; step < PROPAGATION_STEPS; step++) {
        // Choose a random existing fracture to propagate from
        if (m_fractures.empty()) break;
        
        int fracIndex = std::uniform_int_distribution<int>(0, static_cast<int>(m_fractures.size()) - 1)(m_engine);
        Fracture& fracture = m_fractures[fracIndex];
        
        // Decide which end to propagate from
        bool fromStart = std::uniform_real_distribution<float>(0.0f, 1.0f)(m_engine) < 0.5f;
        Math::Vector3 propPoint = fromStart ? fracture.start : fracture.end;
        Math::Vector3 dirVector = fromStart ? (fracture.start - fracture.end).normalized() : 
                                            (fracture.end - fracture.start).normalized();
        
        // Calculate stress field direction (perpendicular to current direction with some randomness)
        Math::Vector3 stressDir = propPoint.cross(dirVector).normalized();
        
        // Add randomization based on irregularity parameter
        stressDir = (stressDir + Math::GeoMath::randomUnitVector(m_engine) * m_params.irregularity).normalized();
        
        // New propagation direction (follows stress field with some influence from original direction)
        Math::Vector3 newDir = (stressDir * 0.8f + dirVector * 0.2f).normalized();
        
        // Length of new segment varies based on step (shorter as we progress)
        float newLength = std::uniform_real_distribution<float>(
            0.1f, 0.3f * (1.0f - static_cast<float>(step) / PROPAGATION_STEPS))(m_engine);
            
        // Calculate new endpoint
        Math::Vector3 newEndPoint = (propPoint.normalized() * std::cos(newLength) + 
                                   newDir * std::sin(newLength)).normalized() * planetRadius;
        
        // Check for intersection with existing fractures
        bool intersects = false;
        for (const auto& f : m_fractures) {
            if (lineSegmentsIntersect(propPoint, newEndPoint, f.start, f.end)) {
                intersects = true;
                break;
            }
        }
        
        // Only add if no intersection found
        if (!intersects) {
            // Width and depth change slightly from parent fracture
            float width = fracture.width * std::uniform_real_distribution<float>(0.8f, 1.2f)(m_engine);
            float depth = fracture.depth * std::uniform_real_distribution<float>(0.7f, 1.1f)(m_engine);
            
            // Add new fracture segment
            m_fractures.push_back({propPoint, newEndPoint, width, depth});
            
            // Sometimes branch in another direction too (creates realistic fracture networks)
            if (std::uniform_real_distribution<float>(0.0f, 1.0f)(m_engine) < 0.3f) {
                // Branch direction is more perpendicular to propagation direction
                Math::Vector3 branchDir = propPoint.cross(newDir).normalized();
                branchDir = (branchDir + Math::GeoMath::randomUnitVector(m_engine) * 0.3f).normalized();
                
                // Branch length is typically shorter
                float branchLength = newLength * std::uniform_real_distribution<float>(0.5f, 0.8f)(m_engine);
                
                // Calculate branch endpoint
                Math::Vector3 branchEndPoint = (propPoint.normalized() * std::cos(branchLength) + 
                                              branchDir * std::sin(branchLength)).normalized() * planetRadius;
                
                // Add branch fracture with slightly different properties
                m_fractures.push_back({
                    propPoint, 
                    branchEndPoint, 
                    width * std::uniform_real_distribution<float>(0.6f, 0.9f)(m_engine),
                    depth * std::uniform_real_distribution<float>(0.7f, 0.9f)(m_engine)
                });
            }
        }
    }
    
    // Apply final irregularity pass
    applyFractalPerturbations();
}

void PlateGenerator::growPlatesFromFractures(int plateCount) {
    // Use the fracture network to define plate boundaries
    // This is a simplified version - in a real implementation, you'd use
    // a more sophisticated algorithm to grow regions
    
    // First, choose seed points for each plate that are away from fractures
    std::vector<Math::Vector3> seedPoints;
    
    const int MAX_ATTEMPTS = 1000;
    
    while (seedPoints.size() < plateCount && seedPoints.size() < MAX_ATTEMPTS) {
        Math::Vector3 candidate = Math::GeoMath::randomUnitVector(m_engine) * m_planetRadius;
        
        // Check if candidate is far enough from all fractures
        bool farEnough = true;
        for (const auto& fracture : m_fractures) {
            float distToFracture = distancePointToLineSegment(
                candidate, fracture.start, fracture.end);
            
            if (distToFracture < fracture.width * m_planetRadius * 2.0f) {
                farEnough = false;
                break;
            }
        }
        
        if (farEnough) {
            // Make sure it's also not too close to other seeds
            for (const auto& existing : seedPoints) {
                if ((existing - candidate).length() < m_planetRadius * 0.3f) {
                    farEnough = false;
                    break;
                }
            }
        }
        
        if (farEnough) {
            seedPoints.push_back(candidate);
        }
    }
    
    // Create plates from seed points
    for (int i = 0; i < seedPoints.size(); i++) {
        Plate plate;
        plate.id = i;
        plate.center = seedPoints[i];
        m_plates.push_back(plate);
    }
    
    // Now create boundaries for each plate
    // In a real implementation, this would be a complex nearest-fracture algorithm
    // This simplified version creates rough boundaries
    const int BOUNDARY_POINTS = 200;
    
    for (auto& plate : m_plates) {
        // Sample points on sphere
        for (int i = 0; i < BOUNDARY_POINTS; i++) {
            // Create a boundary point
            float theta = std::uniform_real_distribution<float>(0.0f, static_cast<float>(M_PI))(m_engine);
            float phi = std::uniform_real_distribution<float>(0.0f, static_cast<float>(2*M_PI))(m_engine);
            
            Math::SphericalCoord boundaryPoint = {
                m_planetRadius,
                theta,
                phi
            };
            
            // Convert to cartesian to check which plate this belongs to
            Math::Vector3 cartPoint(
                boundaryPoint.radius * std::sin(boundaryPoint.theta) * std::cos(boundaryPoint.phi),
                boundaryPoint.radius * std::sin(boundaryPoint.theta) * std::sin(boundaryPoint.phi),
                boundaryPoint.radius * std::cos(boundaryPoint.theta)
            );
            
            // Find the closest seed point, but also consider fractures
            int closestPlateId = -1;
            float minDistance = std::numeric_limits<float>::max();
            
            for (const auto& otherPlate : m_plates) {
                float dist = (cartPoint - otherPlate.center).length();
                
                // Check if there's a fracture between this point and the plate center
                bool fractureInBetween = false;
                for (const auto& fracture : m_fractures) {
                    // Simplified check - in reality this is more complex
                    if (lineSegmentsIntersect(
                            cartPoint, otherPlate.center,
                            fracture.start, fracture.end)) {
                        fractureInBetween = true;
                        break;
                    }
                }
                
                // If there's no fracture in between and it's closer
                if (!fractureInBetween && dist < minDistance) {
                    minDistance = dist;
                    closestPlateId = otherPlate.id;
                }
            }
            
            // Add this point to the correct plate
            if (closestPlateId == plate.id) {
                plate.boundaries.push_back(boundaryPoint);
            }
        }
    }
}

float PlateGenerator::getRandomContinentalThickness() {
    return std::uniform_real_distribution<float>(70.0f, 100.0f)(m_engine);
}

float PlateGenerator::getRandomOceanicThickness() {
    return std::uniform_real_distribution<float>(30.0f, 60.0f)(m_engine);
}

float PlateGenerator::getRandomContinentalDensity() {
    return std::uniform_real_distribution<float>(2700.0f, 2900.0f)(m_engine);
}

float PlateGenerator::getRandomOceanicDensity() {
    return std::uniform_real_distribution<float>(3000.0f, 3300.0f)(m_engine);
}

void PlateGenerator::assignPlateProperties(bool pangea) {
    // Assign realistic properties to plates
    
    // Calculate cluster center for Pangea (if enabled)
    Math::Vector3 pangeaCenter;
    if (pangea) {
        pangeaCenter = Math::GeoMath::randomUnitVector(m_engine) * m_planetRadius;
    }
    
    for (auto& plate : m_plates) {
        // First decide if continental or oceanic
        // For Pangea, plates closer to the center are more likely to be continental
        bool isContinental;
        
        if (pangea) {
            float distToPangeaCenter = (plate.center - pangeaCenter).length() / m_planetRadius;
            float continentalProbability = std::max(0.0f, 1.0f - distToPangeaCenter);
            isContinental = std::uniform_real_distribution<float>(0.0f, 1.0f)(m_engine) < continentalProbability;
        } else {
            // Random distribution across the planet
            isContinental = std::uniform_real_distribution<float>(0.0f, 1.0f)(m_engine) < m_params.continentalPercent;
        }
        
        // Set properties based on crust type
        plate.isOceanic = !isContinental;
        
        if (isContinental) {
            // Continental crust: thicker, less dense
            plate.thickness = getRandomContinentalThickness();
            plate.density = getRandomContinentalDensity();
            plate.buoyancy = 1.2f;
            plate.age = std::uniform_real_distribution<float>(50.0f, 3000.0f)(m_engine); // Continental crust can be very old
        } else {
            // Oceanic crust: thinner, denser
            plate.thickness = getRandomOceanicThickness();
            plate.density = getRandomOceanicDensity();
            plate.buoyancy = 0.8f;
            plate.age = std::uniform_real_distribution<float>(0.0f, 180.0f)(m_engine); // Oceanic crust is younger
        }
        
        // Calculate velocity based on position and crust type
        plate.velocity = calculateInitialVelocity(plate.center, isContinental);
    }
    
    // Identify neighboring plates (needed for boundary interactions)
    for (auto& plate : m_plates) {
        std::unordered_set<int> neighbors;
        
        // For each boundary point, check if it's close to boundaries of other plates
        for (const auto& boundary : plate.boundaries) {
            Math::Vector3 boundaryPoint(
                boundary.radius * std::sin(boundary.theta) * std::cos(boundary.phi),
                boundary.radius * std::sin(boundary.theta) * std::sin(boundary.phi),
                boundary.radius * std::cos(boundary.theta)
            );
            
            for (const auto& otherPlate : m_plates) {
                if (otherPlate.id == plate.id) continue;
                
                for (const auto& otherBoundary : otherPlate.boundaries) {
                    Math::Vector3 otherPoint(
                        otherBoundary.radius * std::sin(otherBoundary.theta) * std::cos(otherBoundary.phi),
                        otherBoundary.radius * std::sin(otherBoundary.theta) * std::sin(otherBoundary.phi),
                        otherBoundary.radius * std::cos(otherBoundary.theta)
                    );
                    
                    if ((boundaryPoint - otherPoint).length() < m_planetRadius * 0.1f) {
                        neighbors.insert(otherPlate.id);
                        break;
                    }
                }
            }
        }
        
        // Store neighbor IDs
        plate.neighborIds.assign(neighbors.begin(), neighbors.end());
    }
}

void PlateGenerator::createWeakZones() {
    // Create weak zones where plates might split in the future
    // This is a simplified version - real implementation would be more complex
    
    // For each plate, potentially create a weak zone
    for (auto& plate : m_plates) {
        // Only some plates get weak zones
        if (std::uniform_real_distribution<float>(0.0f, 1.0f)(m_engine) < 0.3f) {
            // Create a fracture through the plate
            Math::Vector3 start = Math::GeoMath::randomUnitVector(m_engine) * m_planetRadius;
            Math::Vector3 end = Math::GeoMath::randomUnitVector(m_engine) * m_planetRadius;
            
            // Random width and depth
            float width = std::uniform_real_distribution<float>(0.005f, 0.02f)(m_engine);
            float depth = std::uniform_real_distribution<float>(5.0f, 20.0f)(m_engine);
            
            // Add fracture
            m_fractures.push_back({start, end, width, depth});
        }
    }
}

void PlateGenerator::balanceCrustTypes(float continentalPercent) {
    // Ensure that continental crust covers approximately the target percentage of the surface
    
    // Count current distribution
    int continentalCount = 0;
    for (const auto& plate : m_plates) {
        if (!plate.isOceanic) continentalCount++;
    }
    
    // Calculate target count
    int targetContinental = static_cast<int>(m_plates.size() * continentalPercent);
    
    // Adjust if needed
    if (continentalCount > targetContinental) {
        // Too many continental plates, convert some to oceanic
        int toConvert = continentalCount - targetContinental;
        
        // Sort plates by size (number of boundary points as proxy)
        std::vector<std::pair<int, size_t>> platesBySize;
        for (int i = 0; i < m_plates.size(); i++) {
            if (!m_plates[i].isOceanic) {
                platesBySize.push_back({i, m_plates[i].boundaries.size()});
            }
        }
        
        // Sort by ascending size
        std::sort(platesBySize.begin(), platesBySize.end(), 
                 [](const auto& a, const auto& b) { return a.second < b.second; });
        
        // Convert smallest continental plates to oceanic
        for (int i = 0; i < std::min(toConvert, static_cast<int>(platesBySize.size())); i++) {
            int plateIdx = platesBySize[i].first;
            m_plates[plateIdx].isOceanic = true;
            m_plates[plateIdx].thickness = getRandomOceanicThickness();
            m_plates[plateIdx].density = getRandomOceanicDensity();
            m_plates[plateIdx].buoyancy = 0.8f;
            m_plates[plateIdx].age = std::uniform_real_distribution<float>(0.0f, 180.0f)(m_engine);
        }
    } 
    else if (continentalCount < targetContinental) {
        // Too few continental plates, convert some oceanic to continental
        int toConvert = targetContinental - continentalCount;
        
        // Sort oceanic plates by size
        std::vector<std::pair<int, size_t>> platesBySize;
        for (int i = 0; i < m_plates.size(); i++) {
            if (m_plates[i].isOceanic) {
                platesBySize.push_back({i, m_plates[i].boundaries.size()});
            }
        }
        
        // Sort by descending size
        std::sort(platesBySize.begin(), platesBySize.end(), 
                 [](const auto& a, const auto& b) { return a.second > b.second; });
        
        // Convert largest oceanic plates to continental
        for (int i = 0; i < std::min(toConvert, static_cast<int>(platesBySize.size())); i++) {
            int plateIdx = platesBySize[i].first;
            m_plates[plateIdx].isOceanic = false;
            m_plates[plateIdx].thickness = getRandomContinentalThickness();
            m_plates[plateIdx].density = getRandomContinentalDensity();
            m_plates[plateIdx].buoyancy = 1.2f;
            m_plates[plateIdx].age = std::uniform_real_distribution<float>(50.0f, 3000.0f)(m_engine);
        }
    }
}

void PlateGenerator::createHotspots() {
    // Create mantle hotspots/plumes
    int numHotspots = static_cast<int>(4 * M_PI * m_planetRadius * m_planetRadius * m_params.hotspotDensity);
    
    for (int i = 0; i < numHotspots; i++) {
        m_hotspots.push_back(Math::GeoMath::randomUnitVector(m_engine) * m_planetRadius);
    }
}

Math::Vector3 PlateGenerator::calculateInitialVelocity(const Math::Vector3& position, bool isContinental) {
    // Create Earth-like velocity patterns
    // Plates move differently based on their position (latitude)
    
    // Extract latitude information from position (y-axis in spherical coords)
    float latitude = std::asin(position.normalized().z);
    
    // Direction changes based on latitude (mimicking Earth's major currents)
    Math::Vector3 direction;
    
    // Simplified Earth plate motion: 
    // - Equatorial regions tend to move east-west
    // - Mid latitudes tend to have more complex movements
    // - Polar regions have different movement patterns
    if (std::abs(latitude) < 0.2f) { // Equatorial
        direction = Math::Vector3(
            std::sin(position.y * 2.0f), // East-west oscillation
            -std::cos(position.x * 2.0f), // North-south component
            0.0f  // Minimal vertical
        ).normalized();
    } 
    else if (std::abs(latitude) < 0.6f) { // Mid-latitudes
        if (latitude > 0) { // Northern hemisphere
            direction = Math::Vector3(
                std::cos(position.y * 3.0f),
                std::sin(position.x * 3.0f),
                0.1f * std::sin(position.x + position.y)
            ).normalized();
        } else { // Southern hemisphere
            direction = Math::Vector3(
                -std::cos(position.y * 3.0f),
                -std::sin(position.x * 3.0f),
                0.1f * std::sin(position.x + position.y)
            ).normalized();
        }
    } 
    else { // Polar regions
        direction = Math::Vector3(
            std::cos(position.x * 4.0f),
            std::sin(position.y * 4.0f),
            0.05f * std::cos(position.x * position.y)
        ).normalized();
    }
    
    // Calculate speed based on plate type
    // Continental plates generally move slower than oceanic plates
    float speed;
    if (isContinental) {
        speed = std::uniform_real_distribution<float>(0.2f, 0.5f)(m_engine);
    } else {
        speed = std::uniform_real_distribution<float>(0.4f, 0.8f)(m_engine);
    }
    
    return direction * speed;
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

// Helper methods for fracture calculations
float PlateGenerator::distancePointToLineSegment(
    const Math::Vector3& point, 
    const Math::Vector3& lineStart, 
    const Math::Vector3& lineEnd) {
    
    Math::Vector3 ab = lineEnd - lineStart;
    Math::Vector3 ac = point - lineStart;
    
    float proj = ac.dot(ab) / ab.lengthSquared();
    
    if (proj < 0) {
        return (point - lineStart).length();
    } else if (proj > 1) {
        return (point - lineEnd).length();
    } else {
        Math::Vector3 closest = lineStart + ab * proj;
        return (point - closest).length();
    }
}

bool PlateGenerator::lineSegmentsIntersect(
    const Math::Vector3& p1, const Math::Vector3& p2,
    const Math::Vector3& p3, const Math::Vector3& p4) {
    
    // This is a simplified approximation for 3D line segment intersection
    // Accurate 3D line segment intersection is more complex
    
    // Create planes containing each line segment and the origin
    Math::Vector3 n1 = p1.cross(p2).normalized();
    Math::Vector3 n2 = p3.cross(p4).normalized();
    
    // Check if the planes are nearly parallel
    if (std::abs(n1.dot(n2)) > 0.9999f) {
        return false;
    }
    
    // Check if p3 and p4 are on opposite sides of plane 1
    float d3 = p3.dot(n1);
    float d4 = p4.dot(n1);
    if (d3 * d4 > 0) {
        return false;
    }
    
    // Check if p1 and p2 are on opposite sides of plane 2
    float d1 = p1.dot(n2);
    float d2 = p2.dot(n2);
    if (d1 * d2 > 0) {
        return false;
    }
    
    // Lines intersect (simplified approximation)
    return true;
}

const std::vector<Plate>& PlateGenerator::getPlates() const {
    return m_plates;
}

std::vector<Simulation::PlateTectonics::Plate> PlateGenerator::convertToTectonicPlates() const {
    std::vector<Simulation::PlateTectonics::Plate> tectonicPlates;
    for(const auto& plate : m_plates) {
        Simulation::PlateTectonics::Plate tectonicPlate;
        tectonicPlate.id = plate.id;
        
        // Convert the boundaries to float vector
        for(const auto& boundary : plate.boundaries) {
            tectonicPlate.boundaries.push_back(boundary.radius);
            tectonicPlate.boundaries.push_back(boundary.theta);
            tectonicPlate.boundaries.push_back(boundary.phi);
        }
        
        // Calculate angular velocity based on plate movement characteristics
        tectonicPlate.angularVelocity = plate.velocity.length() / (m_planetRadius * 1000.0f); // Convert km to meters
        tectonicPlate.thickness = plate.thickness;
        tectonicPlate.density = plate.density;
        
        // Add new properties
        tectonicPlate.isOceanic = plate.isOceanic;
        tectonicPlate.age = plate.age;
        
        // Add neighbors
        tectonicPlate.neighborIds = plate.neighborIds;
        
        tectonicPlates.push_back(tectonicPlate);
    }
    return tectonicPlates;
}

// Add this to PlateGenerator.cpp
void PlateGenerator::applyFractalPerturbations() {
    // Apply fractal noise to make fractures and plate boundaries more irregular
    for (auto& fracture : m_fractures) {
        // Midpoint displacement algorithm to add irregularity
        
        // Start with original endpoints
        Math::Vector3 start = fracture.start;
        Math::Vector3 end = fracture.end;
        
        // Apply small random perturbations
        fracture.start += Math::GeoMath::randomUnitVector(m_engine) * 
                         (m_params.irregularity * m_planetRadius * 0.02f);
        fracture.end += Math::GeoMath::randomUnitVector(m_engine) * 
                       (m_params.irregularity * m_planetRadius * 0.02f);
        
        // Normalize back to sphere surface
        fracture.start = fracture.start.normalized() * m_planetRadius;
        fracture.end = fracture.end.normalized() * m_planetRadius;
    }
    
    // Also add some irregular perturbations to plate boundaries
    for (auto& plate : m_plates) {
        for (auto& boundary : plate.boundaries) {
            // Add small random perturbations to theta and phi
            float thetaPerturb = std::uniform_real_distribution<float>(
                -0.05f * m_params.irregularity, 
                0.05f * m_params.irregularity)(m_engine);
                
            float phiPerturb = std::uniform_real_distribution<float>(
                -0.05f * m_params.irregularity, 
                0.05f * m_params.irregularity)(m_engine);
                
            boundary.theta += thetaPerturb;
            boundary.phi += phiPerturb;
            
            // Keep theta within valid range [0, π]
            boundary.theta = std::max(0.01f, std::min(static_cast<float>(M_PI) - 0.01f, boundary.theta));
        }
    }
}

} // namespace Tools
} // namespace AeonTerra