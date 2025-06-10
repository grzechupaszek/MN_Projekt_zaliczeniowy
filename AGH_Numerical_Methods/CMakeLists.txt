cmake_minimum_required(VERSION 3.16)
project(AGH_Numerical_Methods
    VERSION 1.0.0
    DESCRIPTION "Biblioteka metod numerycznych dla inżynierii obliczeniowej AGH"
    LANGUAGES CXX
)

# Konfiguracja standardu C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Opcje konfiguracji
option(BUILD_TESTS "Budowanie testów jednostkowych" ON)
option(BUILD_EXAMPLES "Budowanie przykładów użycia" ON)
option(BUILD_DOCUMENTATION "Budowanie dokumentacji" OFF)

# Flagi kompilacji w zależności od trybu
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wextra -Wpedantic")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
endif()

# Ścieżki do katalogów
set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)
set(EXAMPLES_DIR ${CMAKE_CURRENT_SOURCE_DIR}/examples)

# Pliki nagłówkowe biblioteki
set(HEADERS
    ${INCLUDE_DIR}/numerical_methods.h
    ${INCLUDE_DIR}/linear_systems.h
    ${INCLUDE_DIR}/interpolation.h
    ${INCLUDE_DIR}/approximation.h
    ${INCLUDE_DIR}/integration.h
    ${INCLUDE_DIR}/differential_equations.h
    ${INCLUDE_DIR}/nonlinear_equations.h
)

# Pliki źródłowe biblioteki
set(SOURCES
    ${SRC_DIR}/linear_systems.cpp
    ${SRC_DIR}/interpolation.cpp
    ${SRC_DIR}/approximation.cpp
    ${SRC_DIR}/integration.cpp
    ${SRC_DIR}/differential_equations.cpp
    ${SRC_DIR}/nonlinear_equations.cpp
)

# Tworzenie biblioteki statycznej
add_library(agh_numerical_methods STATIC ${SOURCES} ${HEADERS})

# Konfiguracja ścieżek nagłówków
target_include_directories(agh_numerical_methods
    PUBLIC
        $<BUILD_INTERFACE:${INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${SRC_DIR}
)

# Ustawienia kompilatora dla biblioteki
target_compile_features(agh_numerical_methods PUBLIC cxx_std_17)

# Linkowanie z bibliotekami matematycznymi
target_link_libraries(agh_numerical_methods PUBLIC m)

# Kompilacja testów jednostkowych
if(BUILD_TESTS)
    enable_testing()
    
    # Lista testów
    set(TEST_SOURCES
        ${TEST_DIR}/test_linear_systems.cpp
        ${TEST_DIR}/test_interpolation.cpp
        ${TEST_DIR}/test_approximation.cpp
        ${TEST_DIR}/test_integration.cpp
        ${TEST_DIR}/test_differential_equations.cpp
        ${TEST_DIR}/test_nonlinear_equations.cpp
    )
    
    # Tworzenie pliku wykonywalnego dla testów
    add_executable(run_tests
        ${TEST_DIR}/run_tests.cpp
        ${TEST_SOURCES}
    )
    
    target_link_libraries(run_tests PRIVATE agh_numerical_methods)
    target_include_directories(run_tests PRIVATE ${INCLUDE_DIR} ${TEST_DIR})
    
    # Dodanie testów do CTest
    add_test(NAME numerical_methods_tests COMMAND run_tests)
    
    # Ustawienie katalogu roboczego dla testów
    set_tests_properties(numerical_methods_tests PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()

# Kompilacja przykładów użycia
if(BUILD_EXAMPLES)
    # Przykład interpolacji
    add_executable(example_interpolation
        ${EXAMPLES_DIR}/example_interpolation.cpp
    )
    target_link_libraries(example_interpolation PRIVATE agh_numerical_methods)
    target_include_directories(example_interpolation PRIVATE ${INCLUDE_DIR})
    
    # Przykład całkowania numerycznego
    add_executable(example_integration
        ${EXAMPLES_DIR}/example_integration.cpp
    )
    target_link_libraries(example_integration PRIVATE agh_numerical_methods)
    target_include_directories(example_integration PRIVATE ${INCLUDE_DIR})
    
    # Przykład równań różniczkowych
    add_executable(example_differential_equations
        ${EXAMPLES_DIR}/example_differential_equations.cpp
    )
    target_link_libraries(example_differential_equations PRIVATE agh_numerical_methods)
    target_include_directories(example_differential_equations PRIVATE ${INCLUDE_DIR})
endif()

# Instalacja biblioteki
install(TARGETS agh_numerical_methods
    EXPORT agh_numerical_methods_targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

# Instalacja nagłówków
install(FILES ${HEADERS}
    DESTINATION include/agh_numerical_methods
)

# Eksport celów dla innych projektów CMake
install(EXPORT agh_numerical_methods_targets
    FILE agh_numerical_methods_targets.cmake
    NAMESPACE agh_numerical_methods::
    DESTINATION lib/cmake/agh_numerical_methods
)

# Generowanie pliku konfiguracyjnego
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/agh_numerical_methods_config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/agh_numerical_methods_config.cmake"
    INSTALL_DESTINATION lib/cmake/agh_numerical_methods
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/agh_numerical_methods_config_version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/agh_numerical_methods_config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/agh_numerical_methods_config_version.cmake"
    DESTINATION lib/cmake/agh_numerical_methods
)

# Informacje o konfiguracji
message(STATUS "Konfiguracja AGH Numerical Methods Library:")
message(STATUS "  Wersja: ${PROJECT_VERSION}")
message(STATUS "  Standard C++: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Tryb budowania: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Testy jednostkowe: ${BUILD_TESTS}")
message(STATUS "  Przykłady: ${BUILD_EXAMPLES}")
message(STATUS "  Dokumentacja: ${BUILD_DOCUMENTATION}")

# Wyświetlenie ścieżek kompilacji
message(STATUS "  Katalog nagłówków: ${INCLUDE_DIR}")
message(STATUS "  Katalog źródeł: ${SRC_DIR}")
message(STATUS "  Flagi kompilatora: ${CMAKE_CXX_FLAGS}")

# Cel do czyszczenia dodatkowych plików
add_custom_target(clean-all
    COMMAND ${CMAKE_BUILD_TOOL} clean
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/Makefile
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake
    COMMENT "Czyszczenie wszystkich plików budowania"
)

# Cel do uruchomienia wszystkich testów z dodatkową diagnostyką
add_custom_target(test-verbose
    COMMAND ${CMAKE_CTEST_COMMAND} --verbose --output-on-failure
    DEPENDS run_tests
    COMMENT "Uruchomienie testów z pełną diagnostyką"
)