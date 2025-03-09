#include "GeologicalTimeline.h"
#include <algorithm>
#include <random>

namespace AeonTerra {
namespace Simulation {

GeologicalTimeline::GeologicalTimeline() {
    // Default constructor creates an empty timeline
}

void GeologicalTimeline::InitializeEarthLike() {
    m_events.clear();
    CreateEarthEvents();
}

void GeologicalTimeline::InitializeFantasy(uint32_t seed) {
    m_events.clear();
    CreateFantasyEvents(seed);
}

void GeologicalTimeline::AddEvent(float timeMA, GeologicalEvent type, 
                                 const std::string& name, const std::string& description,
                                 const std::function<void()>& callback) {
    GeologicalTimelineEvent event;
    event.timeMA = timeMA;
    event.type = type;
    event.name = name;
    event.description = description;
    event.callback = callback;
    
    m_events.push_back(event);
    
    // Keep events sorted from oldest to newest
    std::sort(m_events.begin(), m_events.end(), 
              [](const GeologicalTimelineEvent& a, const GeologicalTimelineEvent& b) {
                  return a.timeMA > b.timeMA;
              });
}

std::vector<GeologicalTimelineEvent> GeologicalTimeline::GetEventsInRange(float startMA, float endMA) const {
    // Make sure startMA is the older time and endMA is the newer
    if (startMA < endMA) {
        std::swap(startMA, endMA);
    }
    
    std::vector<GeologicalTimelineEvent> result;
    
    for (const auto& event : m_events) {
        if (event.timeMA <= startMA && event.timeMA >= endMA) {
            result.push_back(event);
        }
    }
    
    return result;
}

std::vector<GeologicalTimelineEvent> GeologicalTimeline::AdvanceTime(float currentTimeMA, float newTimeMA) {
    // Make sure currentTimeMA is the older time and newTimeMA is the newer
    if (currentTimeMA < newTimeMA) {
        std::swap(currentTimeMA, newTimeMA);
    }
    
    std::vector<GeologicalTimelineEvent> triggeredEvents;
    
    // Find all events in the time range
    auto events = GetEventsInRange(currentTimeMA, newTimeMA);
    
    // Process each event and its callback
    for (const auto& event : events) {
        triggeredEvents.push_back(event);
        
        // Call the callback if defined
        if (event.callback) {
            event.callback();
        }
    }
    
    return triggeredEvents;
}

std::string GeologicalTimeline::GetGeologicalPeriodName(float timeMA) const {
    // Return the geological period name based on time
    if (timeMA > TIME_EARTH_FORMATION) return "Pre-Formation";
    if (timeMA > TIME_PERMIAN_START) return "Paleozoic Era";
    if (timeMA > TIME_TRIASSIC_START) return "Permian Period";
    if (timeMA > TIME_JURASSIC_START) return "Triassic Period";
    if (timeMA > TIME_CRETACEOUS_START) return "Jurassic Period";
    if (timeMA > TIME_PALEOGENE_START) return "Cretaceous Period";
    if (timeMA > TIME_NEOGENE_START) return "Paleogene Period";
    if (timeMA > TIME_QUATERNARY_START) return "Neogene Period";
    return "Quaternary Period";
}

void GeologicalTimeline::CreateEarthEvents() {
    // Add Earth's major geological events
    AddEvent(4500.0f, GeologicalEvent::PlanetFormation, "Earth Formation", 
             "Formation of planet Earth from accretion of material in solar system");
    
    AddEvent(4450.0f, GeologicalEvent::CoreSolidification, "Core Formation", 
             "Differentiation of Earth's iron core from silicate mantle");
    
    AddEvent(4400.0f, GeologicalEvent::CrustFormation, "Crust Formation", 
             "Formation of Earth's first solid crust");
    
    AddEvent(3800.0f, GeologicalEvent::OceanFormation, "First Oceans", 
             "Formation of Earth's first oceans as planet cooled");
    
    AddEvent(3000.0f, GeologicalEvent::PlateTectonicsStart, "Plate Tectonics Begin", 
             "Start of modern-style plate tectonics");
    
    AddEvent(250.0f, GeologicalEvent::SupercontinentFormation, "Pangea Formation", 
             "Assembly of most recent supercontinent Pangea");
    
    AddEvent(200.0f, GeologicalEvent::SupercontinentBreakup, "Pangea Breakup", 
             "Fragmentation of Pangea supercontinent");
    
    AddEvent(66.0f, GeologicalEvent::MajorExtinction, "Cretaceous-Paleogene Extinction", 
             "Extinction event that eliminated non-avian dinosaurs");
}

void GeologicalTimeline::CreateFantasyEvents(uint32_t seed) {
    // Create a randomized timeline for fantasy worlds
    std::mt19937 rng(seed);
    
    // Common distributions for time periods
    std::uniform_real_distribution<float> planetFormationTime(4000.0f, 6000.0f);
    std::uniform_real_distribution<float> earlyEraJitter(0.0f, 500.0f);
    
    // Create base planetary timeline
    float formationTime = planetFormationTime(rng);
    AddEvent(formationTime, GeologicalEvent::PlanetFormation, "Planet Formation", 
             "Formation of fantasy world from cosmic materials");
    
    AddEvent(formationTime - 50.0f - earlyEraJitter(rng), GeologicalEvent::CoreSolidification, "Core Formation", 
             "Differentiation of planetary core");
    
    AddEvent(formationTime - 100.0f - earlyEraJitter(rng), GeologicalEvent::CrustFormation, "Crust Formation", 
             "Formation of first solid crust");
    
    AddEvent(formationTime - 500.0f - earlyEraJitter(rng), GeologicalEvent::OceanFormation, "First Oceans", 
             "Formation of primordial oceans");
    
    // Randomize number of supercontinents (2-5)
    std::uniform_int_distribution<int> numSupercontinents(2, 5);
    int superContinentCount = numSupercontinents(rng);
    float lastBreakupTime = formationTime - 1000.0f;
    
    for (int i = 0; i < superContinentCount; i++) {
        std::string name = "Ancient Supercontinent " + std::to_string(i+1);
        if (i == superContinentCount - 1) {
            name = "Recent Supercontinent";
        }
        
        // Formation time (older to newer)
        float formationJitter = earlyEraJitter(rng);
        float formationTime = lastBreakupTime - 200.0f - formationJitter;
        
        AddEvent(formationTime, GeologicalEvent::SupercontinentFormation, name + " Formation", 
                 "Assembly of major supercontinent");
        
        // Breakup time
        float breakupJitter = earlyEraJitter(rng);
        float breakupTime = formationTime - 300.0f - breakupJitter;
        
        AddEvent(breakupTime, GeologicalEvent::SupercontinentBreakup, name + " Breakup", 
                 "Fragmentation of supercontinent");
                 
        lastBreakupTime = breakupTime;
    }
}

} // namespace Simulation
} // namespace AeonTerra