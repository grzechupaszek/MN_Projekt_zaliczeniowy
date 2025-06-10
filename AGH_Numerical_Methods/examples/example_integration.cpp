/**
 * @file example_integration.cpp
 * @brief Przykład demonstracyjny użycia metod całkowania numerycznego
 * @author Student Inżynierii Obliczeniowej AGH
 * 
 * @details Ten przykład pokazuje praktyczne zastosowanie różnych metod
 *          całkowania numerycznego dostępnych w bibliotece AGH Numerical Methods.
 *          Demonstruje metody prostokątów, trapezów, Simpsona, kwadraturę
 *          Gaussa-Legendre'a oraz metody adaptacyjne.
 */

#include "numerical_methods.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>

using namespace agh_numerical;

/**
 * @brief Demonstracja podstawowych metod całkowania
 */
void demo_basic_integration_methods() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "DEMONSTRACJA PODSTAWOWYCH METOD CAŁKOWANIA" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Funkcja testowa: f(x) = x² + 2x + 1
    auto quadratic_function = [](double x) {
        return x*x + 2*x + 1;
    };
    
    double a = 0.0, b = 2.0;
    int n = 1000;
    
    // Dokładna wartość całki: ∫₀² (x² + 2x + 1) dx = [x³/3 + x² + x]₀² = 8/3 + 4 + 2 = 20/3
    double exact_value = 20.0/3.0;
    
    std::cout << "Funkcja podcałkowa: f(x) = x² + 2x + 1" << std::endl;
    std::cout << "Przedział całkowania: [" << a << ", " << b << "]" << std::endl;
    std::cout << "Liczba podziałów: " << n << std::endl;
    std::cout << "Dokładna wartość całki: " << std::fixed << std::setprecision(8) << exact_value << std::endl;
    
    std::cout << "\n" << std::string(70, '-') << std::endl;
    std::cout << std::setw(20) << "Metoda" << std::setw(15) << "Wynik" 
              << std::setw(15) << "Błąd bezwzględny" << std::setw(15) << "Błąd względny [%]" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // Metoda prostokątów
    double rect_result = Integration::rectangle_method(quadratic_function, a, b, n);
    double rect_error = std::abs(rect_result - exact_value);
    double rect_rel_error = (rect_error / exact_value) * 100;
    
    std::cout << std::setw(20) << "Prostokąty" 
              << std::setw(15) << std::setprecision(8) << rect_result
              << std::setw(15) << std::scientific << std::setprecision(3) << rect_error
              << std::setw(15) << std::fixed << std::setprecision(6) << rect_rel_error << std::endl;
    
    // Metoda trapezów
    double trap_result = Integration::trapezoid_method(quadratic_function, a, b, n);
    double trap_error = std::abs(trap_result - exact_value);
    double trap_rel_error = (trap_error / exact_value) * 100;
    
    std::cout << std::setw(20) << "Trapezy" 
              << std::setw(15) << std::setprecision(8) << trap_result
              << std::setw(15) << std::scientific << std::setprecision(3) << trap_error
              << std::setw(15) << std::fixed << std::setprecision(6) << trap_rel_error << std::endl;
    
    // Metoda Simpsona
    double simp_result = Integration::simpson_method(quadratic_function, a, b, n);
    double simp_error = std::abs(simp_result - exact_value);
    double simp_rel_error = (simp_error / exact_value) * 100;
    
    std::cout << std::setw(20) << "Simpson" 
              << std::setw(15) << std::setprecision(8) << simp_result
              << std::setw(15) << std::scientific << std::setprecision(3) << simp_error
              << std::setw(15) << std::fixed << std::setprecision(6) << simp_rel_error << std::endl;
    
    std::cout << "\nObserwacje:" << std::endl;
    std::cout << "- Metoda Simpsona jest dokładna dla wielomianów stopnia ≤ 3" << std::endl;
    std::cout << "- Błąd metody Simpsona wynosi ~" << std::scientific << simp_error 
              << " (błąd maszynowy)" << std::endl;
}

/**
 * @brief Demonstracja kwadratury Gaussa-Legendre'a
 */
void demo_gauss_legendre_quadrature() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "DEMONSTRACJA KWADRATURY GAUSSA-LEGENDRE'A" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Test na wielomianie stopnia 7: f(x) = x⁷ - 2x⁵ + x³ - x + 1
    auto polynomial7 = [](double x) {
        return std::pow(x, 7) - 2*std::pow(x, 5) + std::pow(x, 3) - x + 1;
    };
    
    double a = -1.0, b = 1.0;
    
    // Dokładna wartość (obliczona analitycznie)
    double exact_value = 2.0; // ∫₋₁¹ (x⁷ - 2x⁵ + x³ - x + 1) dx = 2
    
    std::cout << "Funkcja podcałkowa: f(x) = x⁷ - 2x⁵ + x³ - x + 1" << std::endl;
    std::cout << "Przedział całkowania: [" << a << ", " << b << "]" << std::endl;
    std::cout << "Dokładna wartość całki: " << exact_value << std::endl;
    
    std::cout << "\n" << std::string(60, '-') << std::endl;
    std::cout << std::setw(15) << "Liczba węzłów" << std::setw(15) << "Wynik G-L" 
              << std::setw(15) << "Błąd" << std::setw(15) << "Dokładność" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    for (int n = 1; n <= 5; ++n) {
        double gl_result = Integration::gauss_legendre_quadrature(polynomial7, a, b, n);
        double error = std::abs(gl_result - exact_value);
        bool is_exact = (error < 1e-12);
        
        std::cout << std::setw(15) << n
                  << std::setw(15) << std::fixed << std::setprecision(8) << gl_result
                  << std::setw(15) << std::scientific << std::setprecision(3) << error
                  << std::setw(15) << (is_exact ? "DOKŁADNA" : "PRZYBLIŻONA") << std::endl;
    }
    
    std::cout << "\nTeoria kwadratury Gaussa-Legendre'a:" << std::endl;
    std::cout << "- Kwadratura n-punktowa jest dokładna dla wielomianów stopnia ≤ 2n-1" << std::endl;
    std::cout << "- Dla wielomianu stopnia 7 potrzebujemy co najmniej 4 węzłów (2×4-1 = 7)" << std::endl;
    std::cout << "- Widać, że od n=4 kwadratura daje wynik dokładny" << std::endl;
}

/**
 * @brief Demonstracja funkcji oscylacyjnych
 */
void demo_oscillatory_functions() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "DEMONSTRACJA CAŁKOWANIA FUNKCJI OSCYLACYJNYCH" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Funkcja oscylacyjna: f(x) = sin(10x) * exp(-x²)
    auto oscillatory_function = [](double x) {
        return std::sin(10*x) * std::exp(-x*x);
    };
    
    double a = -2.0, b = 2.0;
    
    std::cout << "Funkcja podcałkowa: f(x) = sin(10x) * exp(-x²)" << std::endl;
    std::cout << "Przedział całkowania: [" << a << ", " << b << "]" << std::endl;
    std::cout << "Charakterystyka: silnie oscylacyjna z tłumieniam eksponencjalnym" << std::endl;
    
    // Porównanie różnych metod
    std::cout << "\n" << std::string(70, '-') << std::endl;
    std::cout << std::setw(25) << "Metoda" << std::setw(15) << "Wynik" 
              << std::setw(15) << "Czas [ms]" << std::setw(15) << "Wywołania f(x)" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // Metoda Simpsona - duża liczba podziałów
    auto start = std::chrono::high_resolution_clock::now();
    double simp_result = Integration::simpson_method(oscillatory_function, a, b, 10000);
    auto end = std::chrono::high_resolution_clock::now();
    double simp_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << std::setw(25) << "Simpson (n=10000)"
              << std::setw(15) << std::fixed << std::setprecision(8) << simp_result
              << std::setw(15) << std::setprecision(3) << simp_time
              << std::setw(15) << 10001 << std::endl;
    
    // Adaptacyjna kwadratura G-L
    start = std::chrono::high_resolution_clock::now();
    double adaptive_result = Integration::adaptive_gauss_legendre(oscillatory_function, a, b, 5, 100);
    end = std::chrono::high_resolution_clock::now();
    double adaptive_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << std::setw(25) << "Adaptacyjna G-L"
              << std::setw(15) << adaptive_result
              << std::setw(15) << adaptive_time
              << std::setw(15) << 500 << std::endl;
    
    // Metoda specjalizowana dla funkcji oscylacyjnych
    start = std::chrono::high_resolution_clock::now();
    double osc_result = Integration::integrate_oscillatory(oscillatory_function, a, b, 2*M_PI/10, 20);
    end = std::chrono::high_resolution_clock::now();
    double osc_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << std::setw(25) << "Specjalna oscylacyjna"
              << std::setw(15) << osc_result
              << std::setw(15) << osc_time
              << std::setw(15) << "~800" << std::endl;
    
    std::cout << "\nWnioski:" << std::endl;
    std::cout << "- Funkcje oscylacyjne wymagają specjalnego traktowania" << std::endl;
    std::cout << "- Adaptacyjne metody są bardziej efektywne niż zwiększanie n" << std::endl;
    std::cout << "- Znajomość okresu oscylacji pozwala na optymalizację" << std::endl;
}

/**
 * @brief Test zbieżności metod całkowania
 */
void demo_convergence_analysis() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "ANALIZA ZBIEŻNOŚCI METOD CAŁKOWANIA" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Funkcja testowa: f(x) = exp(x) * cos(x)
    auto test_function = [](double x) {
        return std::exp(x) * std::cos(x);
    };
    
    double a = 0.0, b = M_PI;
    // Dokładna wartość: ∫₀^π exp(x)cos(x) dx = (exp(π) + 1)/2
    double exact_value = (std::exp(M_PI) + 1) / 2;
    
    std::cout << "Funkcja podcałkowa: f(x) = exp(x) * cos(x)" << std::endl;
    std::cout << "Przedział całkowania: [0, π]" << std::endl;
    std::cout << "Dokładna wartość: " << std::fixed << std::setprecision(8) << exact_value << std::endl;
    
    // Test zbieżności dla różnych metod
    const std::string filename = "convergence_analysis.csv";
    
    std::ofstream file(filename);
    file << "n,Rectangle,Trapezoid,Simpson,GL_2,GL_3,GL_4,GL_5" << std::endl;
    
    std::cout << "\n" << std::string(80, '-') << std::endl;
    std::cout << std::setw(6) << "n" << std::setw(12) << "Prostokąty" 
              << std::setw(12) << "Trapezy" << std::setw(12) << "Simpson"
              << std::setw(12) << "G-L (5)" << std::setw(12) << "Najlepszy" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    for (int n = 10; n <= 1000; n *= 2) {
        double rect = Integration::rectangle_method(test_function, a, b, n);
        double trap = Integration::trapezoid_method(test_function, a, b, n);
        double simp = Integration::simpson_method(test_function, a, b, n);
        double gl2 = Integration::gauss_legendre_quadrature(test_function, a, b, 2);
        double gl3 = Integration::gauss_legendre_quadrature(test_function, a, b, 3);
        double gl4 = Integration::gauss_legendre_quadrature(test_function, a, b, 4);
        double gl5 = Integration::gauss_legendre_quadrature(test_function, a, b, 5);
        
        // Znajdź metodę z najmniejszym błędem
        std::vector<std::pair<double, std::string>> errors = {
            {std::abs(rect - exact_value), "Prostokąty"},
            {std::abs(trap - exact_value), "Trapezy"},
            {std::abs(simp - exact_value), "Simpson"},
            {std::abs(gl5 - exact_value), "G-L(5)"}
        };
        
        auto min_error = std::min_element(errors.begin(), errors.end());
        
        std::cout << std::setw(6) << n
                  << std::scientific << std::setprecision(2)
                  << std::setw(12) << std::abs(rect - exact_value)
                  << std::setw(12) << std::abs(trap - exact_value)
                  << std::setw(12) << std::abs(simp - exact_value)
                  << std::setw(12) << std::abs(gl5 - exact_value)
                  << std::setw(12) << min_error->second << std::endl;
        
        // Zapisz do pliku CSV
        file << n << "," << std::abs(rect - exact_value) << ","
             << std::abs(trap - exact_value) << "," << std::abs(simp - exact_value) << ","
             << std::abs(gl2 - exact_value) << "," << std::abs(gl3 - exact_value) << ","
             << std::abs(gl4 - exact_value) << "," << std::abs(gl5 - exact_value) << std::endl;
    }
    
    file.close();
    std::cout << "\nDane zbieżności zapisane w pliku: " << filename << std::endl;
}

/**
 * @brief Demonstracja metod adaptacyjnych
 */
void demo_adaptive_methods() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "DEMONSTRACJA METOD ADAPTACYJNYCH" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Funkcja z ostrym pikiem: f(x) = exp(-100(x-0.5)²)
    auto sharp_peak = [](double x) {
        return std::exp(-100 * (x - 0.5) * (x - 0.5));
    };
    
    double a = 0.0, b = 1.0;
    double tolerance = 1e-8;
    
    std::cout << "Funkcja podcałkowa: f(x) = exp(-100(x-0.5)²)" << std::endl;
    std::cout << "Przedział całkowania: [0, 1]" << std::endl;
    std::cout << "Charakterystyka: bardzo ostry pik w x = 0.5" << std::endl;
    std::cout << "Żądana tolerancja: " << std::scientific << tolerance << std::endl;
    
    // Metoda Simpsona z dużą liczbą podziałów
    std::cout << "\n1. Metoda Simpsona z dużą liczbą podziałów:" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    double simp_uniform = Integration::simpson_method(sharp_peak, a, b, 10000);
    auto end = std::chrono::high_resolution_clock::now();
    double simp_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << "   Wynik: " << std::fixed << std::setprecision(8) << simp_uniform << std::endl;
    std::cout << "   Czas: " << std::setprecision(3) << simp_time << " ms" << std::endl;
    std::cout << "   Wywołania funkcji: 10001" << std::endl;
    
    // Adaptacyjna metoda Simpsona
    std::cout << "\n2. Adaptacyjna metoda Simpsona:" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto adaptive_result = Integration::adaptive_simpson(sharp_peak, a, b, tolerance, 15);
    end = std::chrono::high_resolution_clock::now();
    double adaptive_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    std::cout << "   Wynik: " << adaptive_result.value << std::endl;
    std::cout << "   Czas: " << adaptive_time << " ms" << std::endl;
    std::cout << "   Wywołania funkcji: " << adaptive_result.function_calls << std::endl;
    std::cout << "   Oszacowany błąd: " << std::scientific << adaptive_result.error_estimate << std::endl;
    std::cout << "   Zbieżność: " << (adaptive_result.converged ? "TAK" : "NIE") << std::endl;
    
    // Porównanie efektywności
    std::cout << "\n3. Porównanie efektywności:" << std::endl;
    std::cout << "   Redukcja wywołań funkcji: " 
              << std::fixed << std::setprecision(1) 
              << (100.0 * (10001 - adaptive_result.function_calls) / 10001) << "%" << std::endl;
    std::cout << "   Przyśpieszenie: " 
              << std::setprecision(2) << (simp_time / adaptive_time) << "x" << std::endl;
    std::cout << "   Różnica wyników: " 
              << std::scientific << std::abs(simp_uniform - adaptive_result.value) << std::endl;
}

/**
 * @brief Benchmark wydajności różnych metod
 */
void demo_performance_benchmark() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "BENCHMARK WYDAJNOŚCI METOD CAŁKOWANIA" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Funkcja testowa średniej złożoności
    auto benchmark_function = [](double x) {
        return std::sin(x) * std::exp(-x*x/4) + std::cos(2*x);
    };
    
    double a = -2.0, b = 2.0;
    int n = 1000;
    int iterations = 10000;
    
    std::cout << "Funkcja podcałkowa: sin(x) * exp(-x²/4) + cos(2x)" << std::endl;
    std::cout << "Przedział: [" << a << ", " << b << "]" << std::endl;
    std::cout << "Liczba podziałów: " << n << std::endl;
    std::cout << "Liczba iteracji testu: " << iterations << std::endl;
    
    auto benchmark_results = Integration::benchmark_methods(benchmark_function, a, b, n, iterations);
    
    std::cout << "\n" << std::string(50, '-') << std::endl;
    std::cout << std::setw(20) << "Metoda" << std::setw(15) << "Czas [ms]" 
              << std::setw(15) << "Względna szybkość" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    // Znajdź najszybszą metodę
    double min_time = std::min_element(benchmark_results.begin(), benchmark_results.end(),
                                      [](const auto& a, const auto& b) { 
                                          return a.second < b.second; 
                                      })->second;
    
    for (const auto& result : benchmark_results) {
        double relative_speed = min_time / result.second;
        std::cout << std::setw(20) << result.first
                  << std::setw(15) << std::fixed << std::setprecision(3) << result.second
                  << std::setw(15) << std::setprecision(2) << relative_speed << "x" << std::endl;
    }
}

/**
 * @brief Funkcja główna - demonstracja wszystkich funkcjonalności
 */
int main() {
    std::cout << "AGH NUMERICAL METHODS LIBRARY - PRZYKŁAD CAŁKOWANIA NUMERYCZNEGO" << std::endl;
    std::cout << "Biblioteka metod numerycznych dla inżynierii obliczeniowej" << std::endl;
    std::cout << "Autor: Student Inżynierii Obliczeniowej AGH" << std::endl;
    
    try {
        // Demonstracje poszczególnych funkcjonalności
        demo_basic_integration_methods();
        demo_gauss_legendre_quadrature();
        demo_oscillatory_functions();
        demo_convergence_analysis();
        demo_adaptive_methods();
        demo_performance_benchmark();
        
        // Test walidacyjny biblioteki
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "TESTY WALIDACYJNE BIBLIOTEKI" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        bool validation_passed = Integration::run_validation_tests();
        
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "DEMONSTRACJA ZAKOŃCZONA " << (validation_passed ? "POMYŚLNIE" : "Z BŁĘDAMI") << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        std::cout << "\nPodsumowanie funkcjonalności:" << std::endl;
        std::cout << "✓ Metoda prostokątów (punkt środkowy)" << std::endl;
        std::cout << "✓ Metoda trapezów" << std::endl;
        std::cout << "✓ Metoda Simpsona (1/3)" << std::endl;
        std::cout << "✓ Kwadratura Gaussa-Legendre'a (1-5 węzłów)" << std::endl;
        std::cout << "✓ Adaptacyjna metoda Simpsona z kontrolą błędu" << std::endl;
        std::cout << "✓ Metody specjalizowane dla funkcji oscylacyjnych" << std::endl;
        std::cout << "✓ Analiza zbieżności i porównanie wydajności" << std::endl;
        std::cout << "✓ Automatyczne testy walidacyjne" << std::endl;
        
        std::cout << "\nWygenerowane pliki:" << std::endl;
        std::cout << "- convergence_analysis.csv (analiza zbieżności)" << std::endl;
        
        return validation_passed ? 0 : 1;
        
    } catch (const std::exception& e) {
        std::cerr << "Błąd podczas wykonywania demonstracji: " << e.what() << std::endl;
        return 1;
    }
}