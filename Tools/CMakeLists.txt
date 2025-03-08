file(GLOB_RECURSE TOOLS_SOURCES "WorldBuilder/*.cpp" "WorldBuilder/*.h")
add_library(Tools STATIC ${TOOLS_SOURCES})
target_link_libraries(Tools PRIVATE Core)
target_include_directories(Tools PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/WorldBuilder
)