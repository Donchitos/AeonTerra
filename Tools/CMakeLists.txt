﻿file(GLOB_RECURSE TOOLS_SOURCES "WorldBuilder/*.cpp" "WorldBuilder/*.h")
add_library(Tools STATIC ${TOOLS_SOURCES})

# Add STATIC definition to match Core
target_compile_definitions(Tools PRIVATE AEONTERRA_CORE_STATIC)

target_link_libraries(Tools PRIVATE Core)
target_include_directories(Tools PUBLIC
    ${CMAKE_SOURCE_DIR}  # Root directory
    ${CMAKE_CURRENT_SOURCE_DIR}  # Tools directory
)