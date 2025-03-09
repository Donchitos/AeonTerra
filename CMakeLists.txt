cmake_minimum_required(VERSION 3.20)
project(AeonTerra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build as static libraries by default
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Define global export macro for static build
add_compile_definitions(AEONTERRA_CORE_STATIC)

# Add definition for PI constants on MSVC
if(MSVC)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

# Enable OpenMP for parallel processing if available
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Modern target-based approach
add_subdirectory(Core)
add_subdirectory(Tools)

# Add the test application
add_executable(AeonTerraTest Main.cpp)
target_link_libraries(AeonTerraTest PRIVATE Core Tools)

# Add our planet simulator application
add_executable(PlanetSimulator PlanetSimulator.cpp)
target_link_libraries(PlanetSimulator PRIVATE Core Tools)

# If this is a Windows build, set the subsystem to console
if(WIN32)
    set_target_properties(AeonTerraTest PROPERTIES
        WIN32_EXECUTABLE FALSE
    )
    set_target_properties(PlanetSimulator PROPERTIES
        WIN32_EXECUTABLE FALSE
    )
endif()

# Ensure we use C++17 filesystem library
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(PlanetSimulator PRIVATE stdc++fs)
endif()

# Install rules
install(TARGETS AeonTerraTest PlanetSimulator
    RUNTIME DESTINATION bin
)
# Add the new visualizer application
add_executable(PlanetVisualizer PlanetVisualizer.cpp)
target_link_libraries(PlanetVisualizer PRIVATE Core Tools)

if(WIN32)
    set_target_properties(PlanetVisualizer PROPERTIES
        WIN32_EXECUTABLE FALSE
    )
endif()

# Ensure we use C++17 filesystem library
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(PlanetVisualizer PRIVATE stdc++fs)
endif()

# Install the visualizer
install(TARGETS PlanetVisualizer
    RUNTIME DESTINATION bin
)