cmake_minimum_required(VERSION 3.15)
project(HGS_VRPSDPTW)
set(CMAKE_CXX_STANDARD 17)

# List all source files for your solver
set(
        src_files
        Program/Genetic.cpp
        Program/Individual.cpp
        Program/LocalSearch.cpp
        Program/Params.cpp
        Program/Population.cpp
        Program/Split.cpp
        Program/InstanceCVRPLIB.cpp   # <- your modified loader
        Program/AlgorithmParameters.cpp
        Program/C_Interface.cpp
)

if (MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

include_directories(Program)

# Object library for reuse
add_library(objlib OBJECT ${src_files})
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

# Shared lib
add_library(lib SHARED $<TARGET_OBJECTS:objlib>)
set_target_properties(lib PROPERTIES OUTPUT_NAME hgs-vrpspdtw)

# If building standalone (not as subproject)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    # Static lib
    add_library(lib_static STATIC $<TARGET_OBJECTS:objlib>)
    set_target_properties(lib_static PROPERTIES OUTPUT_NAME hgs-vrpspdtw_static)

    # Main executable
    add_executable(bin Program/main.cpp)
    target_link_libraries(bin PRIVATE lib_static)
    set_target_properties(bin PROPERTIES OUTPUT_NAME hgs-vrpspdtw)

    # Install targets
    install(TARGETS lib
            DESTINATION lib)
    install(TARGETS bin
            DESTINATION bin)
    install(FILES Program/AlgorithmParameters.h Program/C_Interface.h
            DESTINATION include)
endif()
