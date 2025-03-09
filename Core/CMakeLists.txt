# Explicitly declare core sources and headers
set(CORE_HEADERS
    CoreExports.h
    Mathematics/Vector3.h
    Mathematics/GeoMath.h
    Simulation/PlateTectonics.h
    Simulation/GeologicalTimeline.h
)

set(CORE_SOURCES
    Simulation/PlateTectonics.cpp
    Simulation/GeologicalTimeline.cpp
)

add_library(Core STATIC ${CORE_HEADERS} ${CORE_SOURCES})

# Configure export macros for Windows DLLs
if(WIN32)
    target_compile_definitions(Core PRIVATE AEONTERRA_CORE_EXPORTS)
    
    # Define M_PI for MSVC which doesn't include it by default
    target_compile_definitions(Core PRIVATE _USE_MATH_DEFINES)
endif()

# Properly expose public headers to dependent projects
target_include_directories(Core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    ${CMAKE_SOURCE_DIR}  # Add root directory to include path
)