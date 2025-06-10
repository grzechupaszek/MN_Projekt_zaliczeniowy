/**
 * @file numerical_methods.h
 * @brief AGH Numerical Methods Library - Główny plik nagłówkowy
 * @author Student Inżynierii Obliczeniowej AGH
 * @version 1.0
 * @date 2025
 * 
 * @details Biblioteka metod numerycznych implementująca algorytmy poznane 
 *          podczas laboratoriów z przedmiotu Metody Numeryczne.
 *          Zawiera implementacje metod dla:
 *          - Rozwiązywania układów równań liniowych
 *          - Interpolacji i aproksymacji
 *          - Całkowania numerycznego
 *          - Rozwiązywania równań różniczkowych
 *          - Rozwiązywania równań nieliniowych
 */

#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include "linear_systems.h"
#include "interpolation.h"
#include "approximation.h"
#include "integration.h"
#include "differential_equations.h"
#include "nonlinear_equations.h"

/**
 * @namespace agh_numerical
 * @brief Przestrzeń nazw dla biblioteki AGH Numerical Methods
 */
namespace agh_numerical {
    
    /**
     * @brief Informacje o bibliotece
     */
    struct LibraryInfo {
        static const char* name() { return "AGH Numerical Methods Library"; }
        static const char* version() { return "1.0.0"; }
        static const char* author() { return "Student Inżynierii Obliczeniowej AGH"; }
        static const char* description() { 
            return "Biblioteka metod numerycznych dla inżynierii obliczeniowej"; 
        }
    };
    
    /**
     * @brief Wypisuje informacje o bibliotece
     */
    void print_library_info();
    
} // namespace agh_numerical

#endif // NUMERICAL_METHODS_H