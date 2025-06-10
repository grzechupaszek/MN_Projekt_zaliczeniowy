/**
 * @file numerical_methods.cpp
 * @brief Implementacja głównych funkcji biblioteki AGH Numerical Methods
 * @author Student Inżynierii Obliczeniowej AGH
 */

#include "numerical_methods.h"
#include <iostream>
#include <iomanip>

namespace agh_numerical {

void print_library_info() {
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "  " << LibraryInfo::name() << std::endl;
    std::cout << "  Wersja: " << LibraryInfo::version() << std::endl;
    std::cout << "  Autor: " << LibraryInfo::author() << std::endl;
    std::cout << "  Opis: " << LibraryInfo::description() << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    std::cout << "\nModuły biblioteki:" << std::endl;
    std::cout << "  • Układy równań liniowych (linear_systems.h)" << std::endl;
    std::cout << "    - Eliminacja Gaussa z częściowym pivotingiem" << std::endl;
    std::cout << "    - Rozkład LU z permutacją wierszy" << std::endl;
    std::cout << "    - Obliczanie wyznacznika macierzy" << std::endl;
    
    std::cout << "  • Interpolacja wielomianowa (interpolation.h)" << std::endl;
    std::cout << "    - Metoda Lagrange'a" << std::endl;
    std::cout << "    - Metoda Newtona z różnicami dzielonymi" << std::endl;
    std::cout << "    - Algorytm Hornera" << std::endl;
    
    std::cout << "  • Aproksymacja funkcji (approximation.h)" << std::endl;
    std::cout << "    - Metoda najmniejszych kwadratów" << std::endl;
    std::cout << "    - Aproksymacja wielomianowa i trygonometryczna" << std::endl;
    std::cout << "    - Aproksymacja robustowa" << std::endl;
    
    std::cout << "  • Całkowanie numeryczne (integration.h)" << std::endl;
    std::cout << "    - Metody Newton-Cotes (prostokąty, trapezy, Simpson)" << std::endl;
    std::cout << "    - Kwadratura Gaussa-Legendre'a" << std::endl;
    std::cout << "    - Metody adaptacyjne z kontrolą błędu" << std::endl;
    
    std::cout << "  • Równania różniczkowe (differential_equations.h)" << std::endl;
    std::cout << "    - Metody jawne (Euler, Heun, RK4)" << std::endl;
    std::cout << "    - Metody adaptacyjne (RK4/5)" << std::endl;
    std::cout << "    - Układy równań i problemy brzegowe" << std::endl;
    
    std::cout << "  • Równania nieliniowe (nonlinear_equations.h)" << std::endl;
    std::cout << "    - Metody zawierające (bisekcja, regula falsi)" << std::endl;
    std::cout << "    - Metody otwarte (Newton, sieczne)" << std::endl;
    std::cout << "    - Optymalizacja i układy równań" << std::endl;
    
    std::cout << "\nCechy biblioteki:" << std::endl;
    std::cout << "  ✓ Wysoka dokładność numeryczna" << std::endl;
    std::cout << "  ✓ Stabilność algorytmów" << std::endl;
    std::cout << "  ✓ Kontrola błędów i walidacja" << std::endl;
    std::cout << "  ✓ Interfejs zgodny z C++17" << std::endl;
    std::cout << "  ✓ Kompleksowe testowanie" << std::endl;
    std::cout << "  ✓ Dokumentacja z przykładami" << std::endl;
    
    std::cout << "\nWspierane kompilatory:" << std::endl;
    std::cout << "  • GCC 7.0+ z obsługą C++17" << std::endl;
    std::cout << "  • Clang 6.0+ z obsługą C++17" << std::endl;
    std::cout << "  • MSVC 2019+ (Visual Studio)" << std::endl;
    
    std::cout << "\nSystemy budowania:" << std::endl;
    std::cout << "  • CMake 3.16+ (zalecany)" << std::endl;
    std::cout << "  • Make z dołączonym Makefile" << std::endl;
    
    std::cout << std::string(60, '-') << std::endl;
    std::cout << "Dla pomocy i dokumentacji odwiedź:" << std::endl;

    std::cout << std::string(60, '=') << std::endl;
}

} // namespace agh_numerical