#pragma once

#include "Core/CoreExports.h"
#include <vector>
#include <string>
#include <functional>

namespace AeonTerra {
namespace Simulation {

// Geological time scale constants (in millions of years ago)
constexpr float TIME_PRESENT_DAY = 0.0f;
constexpr float TIME_QUATERNARY_START = 2.58f;
constexpr float TIME_NEOGENE_START = 23.03f;
constexpr float TIME_PALEOGENE_START = 66.0f;
constexpr float TIME_CRETACEOUS_START = 145.0f;
constexpr float TIME_JURASSIC_START = 201.3f;
constexpr float TIME_TRIASSIC_START = 251.902f;
constexpr float TIME_PERMIAN_START = 298.9f;
constexpr float TIME_EARTH_FORMATION = 4500.0f;

// Key geological events
enum class GeologicalEvent {
    PlanetFormation,           // Formation of planet
    CoreSolidification,        // Core formation and solidification
    CrustFormation,            // Initial crust formation
    OceanFormation,            // First oceans appear
    PlateTectonicsStart,       // Plate tectonics begin
    SupercontinentFormation,   // Formation of a supercontinent (e.g., Pangea)
    SupercontinentBreakup,     // Breakup of a supercontinent
    MajorExtinction,           // Major extinction event
    MountainOrogeny,           // Major mountain building event
    IceAge,                    // Ice age
    ClimateWarming,            // Global warming period
    Custom                     // Custom event
};

// Event structure
struct GeologicalTimelineEvent {
    float timeMA;              // Time in millions of years ago
    GeologicalEvent type;
    std::string name;
    std::string description;
    std::function<void()> callback;  // Function to execute when event occurs
};

// Timeline class
class AEONTERRA_API GeologicalTimeline {
public:
    GeologicalTimeline();
    
    // Initialize with Earth-like timeline
    void InitializeEarthLike();
    
    // Initialize with random fantasy timeline
    void InitializeFantasy(uint32_t seed);
    
    // Add a custom event
    void AddEvent(float timeMA, GeologicalEvent type, 
                 const std::string& name, const std::string& description,
                 const std::function<void()>& callback = nullptr);
    
    // Get events in a time range (inclusive)
    std::vector<GeologicalTimelineEvent> GetEventsInRange(float startMA, float endMA) const;
    
    // Process timeline from current time to next time
    // Returns events triggered during this time range
    std::vector<GeologicalTimelineEvent> AdvanceTime(float currentTimeMA, float newTimeMA);
    
    // Get current geological period name
    std::string GetGeologicalPeriodName(float timeMA) const;
    
private:
    // Events sorted by time (descending, earlier events first)
    std::vector<GeologicalTimelineEvent> m_events;
    
    // Helper to create an Earth-like timeline
    void CreateEarthEvents();
    
    // Helper to create a randomized fantasy timeline
    void CreateFantasyEvents(uint32_t seed);
};

} // namespace Simulation
} // namespace AeonTerra