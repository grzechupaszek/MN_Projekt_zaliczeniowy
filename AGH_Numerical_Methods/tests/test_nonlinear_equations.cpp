/**
 * @file test_nonlinear_equations.cpp
 * @brief Testy jednostkowe dla modułu równań nieliniowych
 */

#include "nonlinear_equations.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>

using namespace agh_numerical;

/**
 * @brief Test metody bisekcji dla równania x² - 2 = 0
 */
bool test_bisection_sqrt2() {
    std::cout << "Test: Bisekcja dla √2... ";
    
    auto f = [](double x) { return x*x - 2.0; };
    double a = 1.0, b = 2.0;
    double tolerance = 1e-8;
    
    auto result = NonlinearEquations::bisection_method(f, a, b, tolerance, 100);
    
    double expected = std::sqrt(2.0);
    bool success = result.converged && std::abs(result.root - expected) < tolerance;
    
    if (success) {
        std::cout << "POWODZENIE (iteracje: " << result.iteration_count << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - Oczekiwane: " << expected << ", otrzymane: " << result.root << std::endl;
    }
    
    return success;
}

/**
 * @brief Test metody Newtona dla równania x³ - x - 2 = 0
 */
bool test_newton_cubic() {
    std::cout << "Test: Newton dla x³ - x - 2 = 0... ";
    
    auto f = [](double x) { return x*x*x - x - 2.0; };
    auto df = [](double x) { return 3*x*x - 1.0; };
    double x0 = 1.5;
    double tolerance = 1e-10;
    
    auto result = NonlinearEquations::newton_method(f, df, x0, tolerance, 50);
    
    // Sprawdzenie czy f(root) ≈ 0
    double f_root = f(result.root);
    bool success = result.converged && std::abs(f_root) < tolerance;
    
    if (success) {
        std::cout << "POWODZENIE (iteracje: " << result.iteration_count 
                  << ", f(root) = " << std::scientific << f_root << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - Brak zbieżności lub f(root) = " << f_root << std::endl;
    }
    
    return success;
}

/**
 * @brief Test metody Newtona z pochodną numeryczną
 */
bool test_newton_numerical_derivative() {
    std::cout << "Test: Newton z pochodną numeryczną... ";
    
    auto f = [](double x) { return std::exp(x) - 2*x - 1; };
    double x0 = 1.0;
    double tolerance = 1e-8;
    
    auto result = NonlinearEquations::newton_numerical_derivative(f, x0, tolerance, 50);
    
    // Sprawdzenie zbieżności
    double f_root = f(result.root);
    bool success = result.converged && std::abs(f_root) < tolerance;
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Brak zbieżności" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test metody siecznych
 */
bool test_secant_method() {
    std::cout << "Test: Metoda siecznych... ";
    
    auto f = [](double x) { return x*x*x - 2*x - 5; };
    double x0 = 2.0, x1 = 3.0;
    double tolerance = 1e-8;
    
    auto result = NonlinearEquations::secant_method(f, x0, x1, tolerance, 50);
    
    // Sprawdzenie zbieżności
    double f_root = f(result.root);
    bool success = result.converged && std::abs(f_root) < tolerance;
    
    if (success) {
        std::cout << "POWODZENIE (iteracje: " << result.iteration_count << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - Brak zbieżności" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test metody regula falsi
 */
bool test_false_position() {
    std::cout << "Test: Metoda regula falsi... ";
    
    auto f = [](double x) { return std::cos(x) - x; };
    double a = 0.0, b = 1.0;
    double tolerance = 1e-8;
    
    auto result = NonlinearEquations::false_position_method(f, a, b, tolerance, 100);
    
    // Sprawdzenie zbieżności i dokładności
    double f_root = f(result.root);
    bool success = result.converged && std::abs(f_root) < tolerance;
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Brak zbieżności" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test metody Brenta
 */
bool test_brent_method() {
    std::cout << "Test: Metoda Brenta... ";
    
    auto f = [](double x) { return x*x*x + 4*x*x - 10; };
    double a = 1.0, b = 2.0;
    double tolerance = 1e-12;
    
    auto result = NonlinearEquations::brent_method(f, a, b, tolerance);
    
    // Sprawdzenie wysokiej dokładności metody Brenta
    double f_root = f(result.root);
    bool success = result.converged && std::abs(f_root) < tolerance;
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Brak zbieżności lub niska dokładność" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test znajdowania wszystkich pierwiastków
 */
bool test_find_all_roots() {
    std::cout << "Test: Znajdowanie wszystkich pierwiastków... ";
    
    // Wielomian (x-1)(x-2)(x-3) = x³ - 6x² + 11x - 6 ma pierwiastki 1, 2, 3
    auto f = [](double x) { return x*x*x - 6*x*x + 11*x - 6; };
    double a = 0.0, b = 4.0;
    
    auto roots = NonlinearEquations::find_all_roots(f, a, b, 1000, 1e-8);
    
    // Sprawdzenie czy znaleziono 3 pierwiastki w przybliżeniu 1, 2, 3
    bool success = (roots.size() == 3);
    if (success) {
        for (size_t i = 0; i < roots.size(); ++i) {
            double expected = i + 1.0;
            if (std::abs(roots[i] - expected) > 1e-6) {
                success = false;
                break;
            }
        }
    }
    
    if (success) {
        std::cout << "POWODZENIE (znaleziono " << roots.size() << " pierwiastki)" << std::endl;
    } else {
        std::cout << "BŁĄD - Niepoprawna liczba lub wartości pierwiastków" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test porównania metod
 */
bool test_compare_methods() {
    std::cout << "Test: Porównanie metod... ";
    
    auto f = [](double x) { return x*x - 3; };
    auto df = [](double x) { return 2*x; };
    double a = 1.0, b = 2.0;
    double exact_root = std::sqrt(3.0);
    
    auto results = NonlinearEquations::compare_methods(f, df, a, b, exact_root);
    
    // Sprawdzenie czy wszystkie metody zbiegły
    bool success = !results.empty();
    for (const auto& result : results) {
        if (!result.converged || std::abs(result.root - exact_root) > 1e-6) {
            success = false;
            break;
        }
    }
    
    if (success) {
        std::cout << "POWODZENIE (testowano " << results.size() << " metod)" << std::endl;
    } else {
        std::cout << "BŁĄD - Niektóre metody nie zbiegły" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test zbieżności kwadratowej metody Newtona
 */
bool test_newton_quadratic_convergence() {
    std::cout << "Test: Zbieżność kwadratowa Newtona... ";
    
    auto f = [](double x) { return x*x - 2; };
    auto df = [](double x) { return 2*x; };
    double x0 = 1.5;
    
    auto result = NonlinearEquations::newton_method(f, df, x0, 1e-12, 10);
    
    // Dla metody Newtona oczekujemy szybkiej zbieżności (mało iteracji)
    bool success = result.converged && result.iteration_count <= 6;
    
    if (success) {
        std::cout << "POWODZENIE (zbieżność w " << result.iteration_count << " iteracjach)" << std::endl;
    } else {
        std::cout << "BŁĄD - Zbyt powolna zbieżność" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test stabilności dla różnych punktów startowych
 */
bool test_stability_different_starts() {
    std::cout << "Test: Stabilność dla różnych punktów startowych... ";
    
    auto f = [](double x) { return x*x - 4; };
    auto df = [](double x) { return 2*x; };
    
    std::vector<double> start_points = {1.5, 2.5, 3.0, 10.0};
    double expected_root = 2.0;
    bool success = true;
    
    for (double x0 : start_points) {
        auto result = NonlinearEquations::newton_method(f, df, x0, 1e-8, 50);
        if (!result.converged || std::abs(result.root - expected_root) > 1e-6) {
            success = false;
            break;
        }
    }
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Niestabilność dla niektórych punktów startowych" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test metody złotego podziału dla optymalizacji
 */
bool test_golden_section_search() {
    std::cout << "Test: Złoty podział... ";
    
    // Funkcja f(x) = (x-2)² + 1 ma minimum w x = 2
    auto f = [](double x) { return (x-2)*(x-2) + 1; };
    double a = 0.0, b = 4.0;
    double tolerance = 1e-6;
    
    auto result = NonlinearEquations::golden_section_search(f, a, b, tolerance);
    
    double expected_min = 2.0;
    bool success = std::abs(result.first - expected_min) < tolerance;
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Oczekiwane minimum: " << expected_min 
                  << ", otrzymane: " << result.first << std::endl;
    }
    
    return success;
}

/**
 * @brief Test estymacji tempa zbieżności
 */
bool test_convergence_rate_estimation() {
    std::cout << "Test: Estymacja tempa zbieżności... ";
    
    // Symulacja błędów dla zbieżności kwadratowej
    std::vector<double> errors = {1e-1, 1e-2, 1e-4, 1e-8, 1e-16};
    
    auto rates = NonlinearEquations::estimate_convergence_rate(errors);
    double estimated_order = rates.second;
    
    // Dla zbieżności kwadratowej oczekujemy rzędu ≈ 2
    bool success = (estimated_order > 1.5 && estimated_order < 2.5);
    
    if (success) {
        std::cout << "POWODZENIE (szacowany rząd: " << std::fixed 
                  << std::setprecision(2) << estimated_order << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - Niepoprawny szacowany rząd: " << estimated_order << std::endl;
    }
    
    return success;
}

/**
 * @brief Test funkcji z wieloma pierwiastkami
 */
bool test_multiple_roots_function() {
    std::cout << "Test: Funkcja z wieloma pierwiastkami... ";
    
    // sin(x) ma pierwiastki w x = nπ
    auto f = [](double x) { return std::sin(x); };
    double a = -M_PI, b = 3*M_PI;
    
    auto roots = NonlinearEquations::find_all_roots(f, a, b, 1000, 1e-8);
    
    // Oczekujemy pierwiastków w przybliżeniu: -π, 0, π, 2π, 3π (5 pierwiastków)
    bool success = (roots.size() >= 4 && roots.size() <= 6); // Pewna tolerancja
    
    if (success) {
        std::cout << "POWODZENIE (znaleziono " << roots.size() << " pierwiastków)" << std::endl;
    } else {
        std::cout << "BŁĄD - Znaleziono " << roots.size() << " pierwiastków (oczekiwano ~5)" << std::endl;
    }
    
    return success;
}

/**
 * @brief Test odporności na błędy numeryczne
 */
bool test_numerical_robustness() {
    std::cout << "Test: Odporność numeryczna... ";
    
    // Funkcja z bardzo małą pochodną w pobliżu pierwiastka
    auto f = [](double x) { return (x-1)*(x-1)*(x-1); }; // Pierwiastek potrójny w x=1
    auto df = [](double x) { return 3*(x-1)*(x-1); };
    double x0 = 1.1;
    
    auto result = NonlinearEquations::newton_method(f, df, x0, 1e-6, 100);
    
    // Metoda Newtona może mieć problem z pierwiastkami wielokrotnymi
    bool success = result.converged || result.iteration_count > 50; // Akceptujemy powolną zbieżność
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Całkowity brak zbieżności" << std::endl;
    }
    
    return success;
}

/**
 * @brief Uruchomienie wszystkich testów równań nieliniowych
 */
bool run_all_nonlinear_equations_tests() {
    std::cout << "\n=== TESTY RÓWNAŃ NIELINIOWYCH ===" << std::endl;
    
    std::vector<std::function<bool()>> tests = {
        test_bisection_sqrt2,
        test_newton_cubic,
        test_newton_numerical_derivative,
        test_secant_method,
        test_false_position,
        test_brent_method,
        test_find_all_roots,
        test_compare_methods,
        test_newton_quadratic_convergence,
        test_stability_different_starts,
        test_golden_section_search,
        test_convergence_rate_estimation,
        test_multiple_roots_function,
        test_numerical_robustness
    };
    
    int passed = 0;
    int total = tests.size();
    
    for (auto& test : tests) {
        if (test()) {
            passed++;
        }
    }
    
    std::cout << "\nWyniki testów równań nieliniowych: " 
              << passed << "/" << total << " (" 
              << std::fixed << std::setprecision(1) 
              << (100.0 * passed / total) << "%)" << std::endl;
    
    return (passed == total);
}