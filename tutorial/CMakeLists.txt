cmake_minimum_required(VERSION 3.15)
project(AriadneConstraintSolvingConstraintSolvingTutorial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(AriadneConstraintSolving REQUIRED)

add_executable(tutorial_constraint_solving tutorial_constraint_solving.cpp)
target_link_libraries(tutorial_constraint_solving PRIVATE AriadneConstraintSolving::ariadne-constraint-solving)
