/**
 * @file example_interpolation.cpp
 * @brief Przykład demonstracyjny użycia metod interpolacji
 * @author Student Inżynierii Obliczeniowej AGH
 * 
 * @details Ten przykład pokazuje praktyczne zastosowanie różnych metod
 *          interpolacji dostępnych w bibliotece AGH Numerical Methods.
 *          Demonstruje interpolację Lagrange'a, Newtona oraz porównanie
 *          wydajności algorytmu Hornera z formą naturalną.
 */

#include "numerical_methods.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>

using namespace agh_numerical;

/**
 * @brief Funkcja demonstrująca podstawowe użycie interpolacji Lagrange'a
 */
void demo_lagrange_interpolation() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "DEMONSTRACJA INTERPOLACJI LAGRANGE'A" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Definicja węzłów interpolacji dla funkcji f(x) = sin(x)
    std::vector<double> nodes = {0.0, M_PI/4, M_PI/2, 3*M_PI/4, M_PI};
    std::vector<double> values;
    
    // Obliczenie wartości funkcji w węzłach
    std::cout << "Węzły interpolacji dla f(x) = sin(x):" << std::endl;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double val = std::sin(nodes[i]);
        values.push_back(val);
        std::cout << "f(" << std::fixed << std::setprecision(4) << nodes[i] 
                  << ") = " << std::setprecision(6) << val << std::endl;
    }
    
    // Interpolacja w różnych punktach
    std::cout << "\nWyniki interpolacji Lagrange'a:" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    std::cout << std::setw(8) << "x" << std::setw(15) << "f(x) dokładne" 
              << std::setw(15) << "f(x) interpolowane" << std::setw(12) << "Błąd" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    for (double x = 0.1; x <= M_PI - 0.1; x += M_PI/10) {
        double exact = std::sin(x);
        double interpolated = Interpolation::lagrange_interpolation(nodes, values, x);
        double error = std::abs(exact - interpolated);
        
        std::cout << std::fixed << std::setprecision(4) << std::setw(8) << x
                  << std::setprecision(8) << std::setw(15) << exact
                  << std::setw(15) << interpolated
                  << std::scientific << std::setw(12) << error << std::endl;
    }
}

/**
 * @brief Demonstracja interpolacji Newtona z różnicami dzielonymi
 */
void demo_newton_interpolation() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "DEMONSTRACJA INTERPOLACJI NEWTONA" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Dane dla funkcji f(x) = x³ - 2x² + x + 1
    std::vector<double> nodes = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> values;
    
    auto cubic_function = [](double x) {
        return x*x*x - 2*x*x + x + 1;
    };
    
    std::cout << "Węzły interpolacji dla f(x) = x³ - 2x² + x + 1:" << std::endl;
    for (size_t i = 0; i < nodes.size(); ++i) {
        double val = cubic_function(nodes[i]);
        values.push_back(val);
        std::cout << "f(" << std::fixed << std::setprecision(1) << nodes[i] 
                  << ") = " << std::setprecision(2) << val << std::endl;
    }
    
    // Obliczenie różnic dzielonych
    auto divided_diffs = Interpolation::compute_divided_differences(nodes, values);
    
    std::cout << "\nTabela różnic dzielonych:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    std::cout << std::setw(8) << "i" << std::setw(8) << "x_i" 
              << std::setw(12) << "f[x_i]" << std::setw(12) << "f[x_i,x_{i+1}]" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        std::cout << std::setw(8) << i << std::setw(8) << nodes[i]
                  << std::setw(12) << std::fixed << std::setprecision(3) << divided_diffs[i][0];
        if (i < nodes.size() - 1) {
            std::cout << std::setw(12) << divided_diffs[i][1];
        }
        std::cout << std::endl;
    }
    
    // Test interpolacji w punkcie x = 2.5
    double test_x = 2.5;
    double exact = cubic_function(test_x);
    double newton_result = Interpolation::newton_interpolation(nodes, divided_diffs, test_x);
    double lagrange_result = Interpolation::lagrange_interpolation(nodes, values, test_x);
    
    std::cout << "\nPorównanie metod dla x = " << test_x << ":" << std::endl;
    std::cout << "Wartość dokładna:      " << std::fixed << std::setprecision(8) << exact << std::endl;
    std::cout << "Interpolacja Newtona:  " << newton_result << std::endl;
    std::cout << "Interpolacja Lagrange: " << lagrange_result << std::endl;
    std::cout << "Błąd Newton:          " << std::scientific << std::abs(exact - newton_result) << std::endl;
    std::cout << "Błąd Lagrange:        " << std::abs(exact - lagrange_result) << std::endl;
}

/**
 * @brief Porównanie wydajności algorytmu Hornera z formą naturalną
 */
void demo_horner_performance() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "PORÓWNANIE WYDAJNOŚCI: HORNER vs FORMA NATURALNA" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Wielomian wysokiego stopnia: P(x) = x^15 + 2x^14 + 3x^13 + ... + 15x + 16
    std::vector<double> coefficients;
    for (int i = 1; i <= 16; ++i) {
        coefficients.push_back(static_cast<double>(i));
    }
    
    std::cout << "Testowany wielomian stopnia " << (coefficients.size() - 1) << std::endl;
    std::cout << "P(x) = 1*x^15 + 2*x^14 + 3*x^13 + ... + 15*x + 16" << std::endl;
    
    double test_x = 1.1;
    int iterations = 1000000;
    
    std::cout << "\nTest wydajności (" << iterations << " iteracji):" << std::endl;
    
    // Pomiar wydajności
    auto timing_results = Interpolation::compare_evaluation_methods(
        coefficients, test_x, iterations);
    
    // Weryfikacja poprawności
    double horner_result = Interpolation::horner_evaluation(coefficients, test_x);
    double natural_result = Interpolation::natural_form_evaluation(coefficients, test_x);
    
    std::cout << "Wyniki obliczeń:" << std::endl;
    std::cout << "  Horner:        " << std::fixed << std::setprecision(12) << horner_result << std::endl;
    std::cout << "  Forma naturalna: " << natural_result << std::endl;
    std::cout << "  Różnica:       " << std::scientific << std::abs(horner_result - natural_result) << std::endl;
    
    std::cout << "\nCzasy wykonania:" << std::endl;
    std::cout << "  Horner:        " << std::fixed << std::setprecision(6) << timing_results.first << " s" << std::endl;
    std::cout << "  Forma naturalna: " << timing_results.second << " s" << std::endl;
    std::cout << "  Przyśpieszenie: " << std::setprecision(2) << (timing_results.second / timing_results.first) << "x" << std::endl;
}

/**
 * @brief Demonstracja optymalizacji liczby węzłów interpolacji
 */
void demo_optimal_nodes_selection() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "OPTYMALIZACJA LICZBY WĘZŁÓW INTERPOLACJI" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Generowanie danych dla funkcji Runge: f(x) = 1/(1 + 25x²)
    std::vector<double> all_nodes, all_values;
    double a = -1.0, b = 1.0;
    int num_points = 50;
    
    auto runge_function = [](double x) {
        return 1.0 / (1.0 + 25*x*x);
    };
    
    std::cout << "Funkcja Runge: f(x) = 1/(1 + 25x²) na przedziale [-1, 1]" << std::endl;
    std::cout << "Generowanie " << num_points << " punktów danych..." << std::endl;
    
    for (int i = 0; i < num_points; ++i) {
        double x = a + i * (b - a) / (num_points - 1);
        all_nodes.push_back(x);
        all_values.push_back(runge_function(x));
    }
    
    // Test różnych liczb węzłów
    std::cout << "\nAnaliza błędu MSE dla różnej liczby węzłów:" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    std::cout << std::setw(12) << "Liczba węzłów" << std::setw(15) << "MSE" << std::setw(15) << "Względny MSE" << std::endl;
    std::cout << std::string(40, '-') << std::endl;
    
    double min_mse = std::numeric_limits<double>::max();
    int optimal_nodes = 5;
    
    for (int num_nodes = 3; num_nodes <= 15; num_nodes += 2) {
        std::vector<double> selected_nodes, selected_values;
        Interpolation::select_nodes(all_nodes, all_values, num_points / num_nodes,
                                   selected_nodes, selected_values);
        
        double mse = Interpolation::compute_mse(all_nodes, all_values,
                                              selected_nodes, selected_values);
        double relative_mse = mse / *std::max_element(all_values.begin(), all_values.end());
        
        std::cout << std::setw(12) << num_nodes
                  << std::scientific << std::setw(15) << mse
                  << std::setw(15) << relative_mse << std::endl;
        
        if (mse < min_mse && mse > 0) {
            min_mse = mse;
            optimal_nodes = num_nodes;
        }
    }
    
    // Znajdowanie optymalnej liczby węzłów automatycznie
    int auto_optimal = Interpolation::find_optimal_nodes_count(all_nodes, all_values, 3, 15);
    
    std::cout << "\nWyniki optymalizacji:" << std::endl;
    std::cout << "  Ręczne znalezienie optimum: " << optimal_nodes << " węzłów (MSE = " 
              << std::scientific << min_mse << ")" << std::endl;
    std::cout << "  Automatyczne znalezienie:   " << auto_optimal << " węzłów" << std::endl;
}

/**
 * @brief Demonstracja wczytywania danych z pliku
 */
void demo_file_input_output() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "DEMONSTRACJA WCZYTYWANIA DANYCH Z PLIKU" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    // Tworzenie przykładowego pliku z danymi
    const std::string filename = "interpolation_example_data.txt";
    
    std::ofstream output_file(filename);
    if (output_file.is_open()) {
        output_file << "xi: 0.0 0.5 1.0 1.5 2.0 2.5 3.0" << std::endl;
        output_file << "f(xi): 1.0 0.8 0.5 0.29 0.14 0.074 0.043" << std::endl;
        output_file.close();
        
        std::cout << "Utworzony plik testowy: " << filename << std::endl;
        std::cout << "Zawartość pliku reprezentuje funkcję f(x) = exp(-x)" << std::endl;
    } else {
        std::cerr << "Nie można utworzyć pliku testowego!" << std::endl;
        return;
    }
    
    // Wczytanie danych z pliku
    std::vector<double> nodes, values;
    bool success = Interpolation::load_data(filename, nodes, values);
    
    if (success) {
        std::cout << "\nDane wczytane pomyślnie:" << std::endl;
        std::cout << "Liczba węzłów: " << nodes.size() << std::endl;
        
        std::cout << "\nWczytane dane:" << std::endl;
        for (size_t i = 0; i < nodes.size(); ++i) {
            std::cout << "  x[" << i << "] = " << std::fixed << std::setprecision(1) 
                      << nodes[i] << ", f(x[" << i << "]) = " 
                      << std::setprecision(3) << values[i] << std::endl;
        }
        
        // Test interpolacji dla wczytanych danych
        double test_point = 1.25;
        double interpolated = Interpolation::lagrange_interpolation(nodes, values, test_point);
        double exact_exp = std::exp(-test_point);
        
        std::cout << "\nTest interpolacji w punkcie x = " << test_point << ":" << std::endl;
        std::cout << "  Interpolowane: " << std::setprecision(6) << interpolated << std::endl;
        std::cout << "  Dokładne exp(-x): " << exact_exp << std::endl;
        std::cout << "  Błąd względny: " << std::setprecision(2) 
                  << std::abs(interpolated - exact_exp) / exact_exp * 100 << "%" << std::endl;
        
        // Generowanie danych do wizualizacji
        std::string plot_filename = "interpolation_plot_data.csv";
        Interpolation::generate_interpolation_data(nodes, values, plot_filename, 100, "lagrange");
        std::cout << "\nDane do wykresu zapisane w pliku: " << plot_filename << std::endl;
        
    } else {
        std::cerr << "Błąd wczytywania danych z pliku!" << std::endl;
    }
}

/**
 * @brief Funkcja główna - demonstracja wszystkich funkcjonalności
 */
int main() {
    std::cout << "AGH NUMERICAL METHODS LIBRARY - PRZYKŁAD INTERPOLACJI" << std::endl;
    std::cout << "Biblioteka metod numerycznych dla inżynierii obliczeniowej" << std::endl;
    std::cout << "Autor: Student Inżynierii Obliczeniowej AGH" << std::endl;
    
    try {
        // Demonstracje poszczególnych funkcjonalności
        demo_lagrange_interpolation();
        demo_newton_interpolation();
        demo_horner_performance();
        demo_optimal_nodes_selection();
        demo_file_input_output();
        
        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "DEMONSTRACJA ZAKOŃCZONA POMYŚLNIE" << std::endl;
        std::cout << std::string(60, '=') << std::endl;
        
        std::cout << "\nPodsumowanie funkcjonalności:" << std::endl;
        std::cout << "✓ Interpolacja Lagrange'a" << std::endl;
        std::cout << "✓ Interpolacja Newtona z różnicami dzielonymi" << std::endl;
        std::cout << "✓ Algorytm Hornera - efektywne obliczanie wielomianów" << std::endl;
        std::cout << "✓ Optymalizacja wyboru węzłów interpolacji" << std::endl;
        std::cout << "✓ Wczytywanie/zapisywanie danych z/do plików" << std::endl;
        std::cout << "✓ Analiza błędów i porównanie metod" << std::endl;
        
        std::cout << "\nWygenerowane pliki:" << std::endl;
        std::cout << "- interpolation_example_data.txt (dane wejściowe)" << std::endl;
        std::cout << "- interpolation_plot_data.csv (dane do wizualizacji)" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Błąd podczas wykonywania demonstracji: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}