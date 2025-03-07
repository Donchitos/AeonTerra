cmake_minimum_required(VERSION 3.21)
project(AeonTerra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(SimCore STATIC
    Core/Simulation/PlateTectonics.cpp
    Core/Mathematics/GeoMath.cpp
    Tools/WorldBuilder/PlateGenerator.cpp
)

target_include_directories(SimCore PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/Core/Simulation
    ${CMAKE_CURRENT_SOURCE_DIR}/Core/Mathematics
    ${CMAKE_CURRENT_SOURCE_DIR}/Tools/WorldBuilder
)