#include "PlateTectonics.h"
#include "../../Tools/WorldBuilder/PlateGenerator.h"
#include <cmath>
#include <algorithm>

namespace AeonTerra {
namespace Simulation {

void PlateTectonics::ResolvePlateCollisions() {
    // Basic collision types: convergent, divergent, transform
    for (auto& plate1 : plates) {
        for (auto& plate2 : plates) {
            if (plate1.id >= plate2.id) continue;

            // Simple density-based collision resolution
            if (plate1.density > plate2.density) {
                // Subduction - plate2 slides under plate1
                plate1.thickness += plate2.thickness * 0.1f;
                plate2.angularVelocity *= 0.8f;
            } else {
                // Mountain building - uplift both plates
                plate1.thickness *= 1.05f;
                plate2.thickness *= 1.05f;
            }
        }
    }
}

} // namespace Simulation
} // namespace AeonTerra