cmake_minimum_required(VERSION 3.15)
project(wing)

find_package(IBAMR 0.18.0 REQUIRED)

set(SRC_FILES
    example.cpp)

add_executable(main3d ${SRC_FILES})

target_link_libraries(main3d 
    PRIVATE 
    IBAMR::IBAMR3d
)

target_compile_definitions(main3d PUBLIC -DNDIM=3)

target_compile_features(main3d PRIVATE cxx_std_11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IBAMR_CXX_FLAGS} -O3 -march=native")