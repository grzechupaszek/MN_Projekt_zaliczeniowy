/**
 * @file test_interpolation.cpp
 * @brief Testy jednostkowe dla modułu interpolacji
 */

#include "interpolation.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace agh_numerical;

/**
 * @brief Test interpolacji Lagrange'a dla wielomianu liniowego
 */
bool test_lagrange_linear() {
    std::cout << "Test: Interpolacja Lagrange'a - funkcja liniowa... ";
    
    // Punkty dla funkcji f(x) = 2x + 1
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {1.0, 3.0, 5.0};
    
    // Test w punkcie x = 1.5, oczekiwana wartość: 2*1.5 + 1 = 4
    double result = Interpolation::lagrange_interpolation(x, y, 1.5);
    double expected = 4.0;
    
    const double tolerance = 1e-12;
    bool success = (std::abs(result - expected) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Oczekiwane: " << expected << ", otrzymane: " << result << std::endl;
    }
    
    return success;
}

/**
 * @brief Test interpolacji Lagrange'a dla wielomianu kwadratowego
 */
bool test_lagrange_quadratic() {
    std::cout << "Test: Interpolacja Lagrange'a - funkcja kwadratowa... ";
    
    // Punkty dla funkcji f(x) = x²
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 4.0};
    
    // Test w punkcie x = 1.5, oczekiwana wartość: 1.5² = 2.25
    double result = Interpolation::lagrange_interpolation(x, y, 1.5);
    double expected = 2.25;
    
    const double tolerance = 1e-12;
    bool success = (std::abs(result - expected) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Oczekiwane: " << expected << ", otrzymane: " << result << std::endl;
    }
    
    return success;
}

/**
 * @brief Test obliczania różnic dzielonych
 */
bool test_divided_differences() {
    std::cout << "Test: Różnice dzielone... ";
    
    // Dane dla funkcji f(x) = x³
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 1.0, 8.0, 27.0};
    
    auto table = Interpolation::compute_divided_differences(x, y);
    
    // Sprawdzenie wymiarów tabeli
    if (table.size() != 4 || table[0].size() != 4) {
        std::cout << "BŁĄD - Niepoprawne wymiary tabeli różnic dzielonych" << std::endl;
        return false;
    }
    
    // Dla funkcji x³, różnice dzielone trzeciego rzędu powinny być stałe = 1
    double third_difference = table[0][3];
    const double tolerance = 1e-10;
    bool success = (std::abs(third_difference - 1.0) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Trzecia różnica dzielona: " << third_difference 
                  << " (oczekiwane: 1.0)" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test interpolacji Newtona
 */
bool test_newton_interpolation() {
    std::cout << "Test: Interpolacja Newtona... ";
    
    // Dane dla funkcji f(x) = x² + x + 1
    std::vector<double> x = {0.0, 1.0, 2.0};
    std::vector<double> y = {1.0, 3.0, 7.0}; // 0+0+1=1, 1+1+1=3, 4+2+1=7
    
    auto divided_diffs = Interpolation::compute_divided_differences(x, y);
    
    // Test w punkcie x = 1.5
    double result = Interpolation::newton_interpolation(x, divided_diffs, 1.5);
    double expected = 1.5*1.5 + 1.5 + 1.0; // 2.25 + 1.5 + 1 = 4.75
    
    const double tolerance = 1e-12;
    bool success = (std::abs(result - expected) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Oczekiwane: " << expected << ", otrzymane: " << result << std::endl;
    }
    
    return success;
}

/**
 * @brief Test zgodności metod Lagrange'a i Newtona
 */
bool test_lagrange_newton_consistency() {
    std::cout << "Test: Zgodność Lagrange vs Newton... ";
    
    // Dane testowe
    std::vector<double> x = {-1.0, 0.0, 1.0, 2.0};
    std::vector<double> y = {0.0, 1.0, 4.0, 9.0}; // Różne wartości
    
    double test_x = 0.5;
    
    // Interpolacja Lagrange'a
    double lagrange_result = Interpolation::lagrange_interpolation(x, y, test_x);
    
    // Interpolacja Newtona
    auto divided_diffs = Interpolation::compute_divided_differences(x, y);
    double newton_result = Interpolation::newton_interpolation(x, divided_diffs, test_x);
    
    const double tolerance = 1e-12;
    bool success = (std::abs(lagrange_result - newton_result) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Lagrange: " << lagrange_result 
                  << ", Newton: " << newton_result << std::endl;
    }
    
    return success;
}

/**
 * @brief Test algorytmu Hornera
 */
bool test_horner_evaluation() {
    std::cout << "Test: Algorytm Hornera... ";
    
    // Wielomian P(x) = 2x³ + 3x² - x + 5
    std::vector<double> coeffs = {2.0, 3.0, -1.0, 5.0}; // Od najwyższej potęgi
    double x = 2.0;
    
    // Oczekiwana wartość: 2(8) + 3(4) - 2 + 5 = 16 + 12 - 2 + 5 = 31
    double result = Interpolation::horner_evaluation(coeffs, x);
    double expected = 31.0;
    
    const double tolerance = 1e-12;
    bool success = (std::abs(result - expected) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Oczekiwane: " << expected << ", otrzymane: " << result << std::endl;
    }
    
    return success;
}

/**
 * @brief Test porównania metod Hornera i formy naturalnej
 */
bool test_horner_vs_natural() {
    std::cout << "Test: Horner vs forma naturalna... ";
    
    // Wielomian testowy
    std::vector<double> coeffs = {1.0, -2.0, 3.0, -4.0, 5.0}; // x⁴ - 2x³ + 3x² - 4x + 5
    double x = 1.5;
    
    double horner_result = Interpolation::horner_evaluation(coeffs, x);
    double natural_result = Interpolation::natural_form_evaluation(coeffs, x);
    
    const double tolerance = 1e-12;
    bool success = (std::abs(horner_result - natural_result) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Horner: " << horner_result 
                  << ", Naturalna: " << natural_result << std::endl;
    }
    
    return success;
}

/**
 * @brief Test wyboru węzłów interpolacji
 */
bool test_node_selection() {
    std::cout << "Test: Wybór węzłów interpolacji... ";
    
    // Wszystkie węzły
    std::vector<double> all_x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> all_y = {0.0, 1.0, 4.0, 9.0, 16.0, 25.0}; // x²
    
    // Wybierz co drugi węzeł (krok = 2)
    std::vector<double> selected_x, selected_y;
    Interpolation::select_nodes(all_x, all_y, 2, selected_x, selected_y);
    
    // Oczekujemy węzły: (0,0), (2,4), (4,16)
    bool success = (selected_x.size() == 3 && selected_y.size() == 3 &&
                   selected_x[0] == 0.0 && selected_y[0] == 0.0 &&
                   selected_x[1] == 2.0 && selected_y[1] == 4.0 &&
                   selected_x[2] == 4.0 && selected_y[2] == 16.0);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Niepoprawny wybór węzłów" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test obliczania błędu MSE
 */
bool test_mse_computation() {
    std::cout << "Test: Obliczanie MSE... ";
    
    // Wszystkie punkty (funkcja liniowa)
    std::vector<double> all_x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> all_y = {0.0, 2.0, 4.0, 6.0, 8.0}; // f(x) = 2x
    
    // Węzły interpolacji (tylko krańce)
    std::vector<double> interp_x = {0.0, 4.0};
    std::vector<double> interp_y = {0.0, 8.0};
    
    // MSE dla interpolacji liniowej przez punkty skrajne
    double mse = Interpolation::compute_mse(all_x, all_y, interp_x, interp_y);
    
    // Dla funkcji liniowej interpolacja przez krańce powinna być dokładna (MSE ≈ 0)
    const double tolerance = 1e-12;
    bool success = (mse < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE (MSE = " << std::scientific << mse << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - MSE = " << mse << " (oczekiwane ≈ 0)" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test znajdowania optymalnej liczby węzłów
 */
bool test_optimal_nodes_count() {
    std::cout << "Test: Optymalna liczba węzłów... ";
    
    // Dane dla funkcji wielomianowej stopnia 3
    std::vector<double> x, y;
    for (int i = 0; i <= 10; ++i) {
        double xi = i * 0.5;
        x.push_back(xi);
        y.push_back(xi*xi*xi - 2*xi*xi + xi + 1); // Wielomian stopnia 3
    }
    
    int optimal_count = Interpolation::find_optimal_nodes_count(x, y, 3, 8);
    
    // Dla wielomianu stopnia 3 optymalna liczba węzłów powinna wynosić 4
    bool success = (optimal_count == 4);
    
    if (success) {
        std::cout << "POWODZENIE (optymalna liczba: " << optimal_count << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - Znaleziono: " << optimal_count << " (oczekiwane: 4)" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test interpolacji dla funkcji sin(x)
 */
bool test_sine_interpolation() {
    std::cout << "Test: Interpolacja sin(x)... ";
    
    // Węzły dla funkcji sin(x) na [0, π]
    std::vector<double> x = {0.0, M_PI/4, M_PI/2, 3*M_PI/4, M_PI};
    std::vector<double> y;
    
    for (double xi : x) {
        y.push_back(std::sin(xi));
    }
    
    // Test w punkcie x = π/3
    double test_x = M_PI / 3;
    double result = Interpolation::lagrange_interpolation(x, y, test_x);
    double expected = std::sin(test_x); // √3/2 ≈ 0.866
    
    // Dla interpolacji funkcji sin(x) oczekujemy dobrej dokładności
    const double tolerance = 0.01; // Większa tolerancja dla funkcji trygonometrycznych
    bool success = (std::abs(result - expected) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Oczekiwane: " << expected << ", otrzymane: " << result 
                  << " (błąd: " << std::abs(result - expected) << ")" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test wydajności algorytmu Hornera
 */
bool test_horner_performance() {
    std::cout << "Test: Wydajność Hornera... ";
    
    // Wielomian wysokiego stopnia
    std::vector<double> coeffs;
    for (int i = 0; i < 20; ++i) {
        coeffs.push_back(1.0 / (i + 1)); // Współczynniki 1, 1/2, 1/3, ...
    }
    
    double x = 0.5;
    int iterations = 10000;
    
    auto timing_results = Interpolation::compare_evaluation_methods(coeffs, x, iterations);
    
    // Horner powinien być szybszy niż forma naturalna
    bool success = (timing_results.first < timing_results.second);
    
    if (success) {
        double speedup = timing_results.second / timing_results.first;
        std::cout << "POWODZENIE (przyśpieszenie: " << std::fixed 
                  << std::setprecision(2) << speedup << "x)" << std::endl;
    } else {
        std::cout << "BŁĄD - Horner nie jest szybszy od formy naturalnej" << std::endl;
    }
    
    return success;
}

/**
 * @brief Uruchomienie wszystkich testów interpolacji
 */
bool run_all_interpolation_tests() {
    std::cout << "\n=== TESTY INTERPOLACJI ===" << std::endl;
    
    std::vector<std::function<bool()>> tests = {
        test_lagrange_linear,
        test_lagrange_quadratic,
        test_divided_differences,
        test_newton_interpolation,
        test_lagrange_newton_consistency,
        test_horner_evaluation,
        test_horner_vs_natural,
        test_node_selection,
        test_mse_computation,
        test_optimal_nodes_count,
        test_sine_interpolation,
        test_horner_performance
    };
    
    int passed = 0;
    int total = tests.size();
    
    for (auto& test : tests) {
        if (test()) {
            passed++;
        }
    }
    
    std::cout << "\nWyniki testów interpolacji: " 
              << passed << "/" << total << " (" 
              << std::fixed << std::setprecision(1) 
              << (100.0 * passed / total) << "%)" << std::endl;
    
    return (passed == total);
}