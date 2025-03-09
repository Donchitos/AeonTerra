cmake_minimum_required(VERSION 3.20)
project(AeonTerra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build as static libraries by default
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Define global export macro for static build
add_compile_definitions(AEONTERRA_CORE_STATIC)

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

# If this is a Windows build, set the subsystem to console
if(WIN32)
    set_target_properties(AeonTerraTest PROPERTIES
        WIN32_EXECUTABLE FALSE
    )
endif()

# Install rules
install(TARGETS AeonTerraTest
    RUNTIME DESTINATION bin
)