# Explicitly declare core sources and headers
set(CORE_HEADERS
    CoreExports.h
    Mathematics/Vector3.h
    Mathematics/GeoMath.h
    Simulation/PlateTectonics.h
    Simulation/PlateTectonics.cpp
)

set(CORE_SOURCES
    ${CORE_HEADERS}
)

add_library(Core STATIC ${CORE_SOURCES})

# Configure export macros for Windows DLLs
if(WIN32)
    target_compile_definitions(Core PRIVATE AEONTERRA_CORE_EXPORTS)
endif()

# Properly expose public headers to dependent projects
target_include_directories(Core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)