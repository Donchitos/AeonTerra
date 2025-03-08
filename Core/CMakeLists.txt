file(GLOB_RECURSE CORE_SOURCES "*.cpp" "*.h")
add_library(Core STATIC ${CORE_SOURCES})
target_compile_definitions(Core PRIVATE AEONTERRA_CORE_STATIC)
target_include_directories(Core PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Mathematics
)