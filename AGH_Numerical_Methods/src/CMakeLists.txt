# Zbierz wszystkie pliki źródłowe
file(GLOB_RECURSE SOURCES
    "approximation.cpp"
    "differential_equations.cpp"
    "integration.cpp"
    "interpolation.cpp"
    "linear_systems.cpp"
    "nonlinear_equations.cpp"
    "numerical_methods.cpp"
)

# Utwórz bibliotekę
add_library(AGH_Numerical_Methods STATIC ${SOURCES})

# Dodaj katalog include do ścieżek wyszukiwania
target_include_directories(AGH_Numerical_Methods PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

# Opcjonalnie: ustaw właściwości biblioteki
set_target_properties(AGH_Numerical_Methods PROPERTIES
    OUTPUT_NAME "agh_numerical_methods"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
