# Dodaj katalog include
include_directories(${CMAKE_SOURCE_DIR}/include)

# Przykłady jako oddzielne pliki wykonywalne
add_executable(example_differential_equations example_differential_equations.cpp)
target_link_libraries(example_differential_equations AGH_Numerical_Methods)

add_executable(example_integration example_integration.cpp)
target_link_libraries(example_integration AGH_Numerical_Methods)

add_executable(example_interpolation example_interpolation.cpp)
target_link_libraries(example_interpolation AGH_Numerical_Methods)
