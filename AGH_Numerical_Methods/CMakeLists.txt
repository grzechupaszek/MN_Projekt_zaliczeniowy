cmake_minimum_required(VERSION 3.16)
project(AGH_Numerical_Methods VERSION 1.0.0 LANGUAGES CXX)

# Ustawienia standardu C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dodaj podkatalogi
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(tests)
