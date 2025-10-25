cmake_minimum_required(VERSION 3.12)
project(StandaloneTransport LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(geometry STATIC
    geometry.cpp
    material.cpp
    raytracer.cpp)
target_include_directories(geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(moc_solver STATIC
    moc_solver.cpp)
target_link_libraries(moc_solver PUBLIC geometry)

add_library(sn_solver STATIC
    sn_solver.cpp)
target_link_libraries(sn_solver PUBLIC geometry)

add_executable(transport_example
    main.cpp)
target_link_libraries(transport_example
    PRIVATE
        moc_solver
        sn_solver
        geometry)
