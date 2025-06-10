/**
 * @file integration.cpp
 * @brief Implementacja metod całkowania numerycznego
 */

#include "integration.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <limits>

namespace agh_numerical {

const double Integration::EPSILON = 1e-15;

double Integration::rectangle_method(std::function<double(double)> f,
                                   double a, double b, int n) {
    if (n <= 0 || a >= b) {
        std::cerr << "Błąd: Niepoprawne parametry całkowania" << std::endl;
        return 0.0;
    }
    
    double h = (b - a) / n;
    double sum = 0.0;
    
    for (int i = 0; i < n; ++i) {
        double x = a + (i + 0.5) * h; // Punkt środkowy
        sum += f(x);
    }
    
    return h * sum;
}

double Integration::trapezoid_method(std::function<double(double)> f,
                                   double a, double b, int n) {
    if (n <= 0 || a >= b) {
        std::cerr << "Błąd: Niepoprawne parametry całkowania" << std::endl;
        return 0.0;
    }
    
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }
    
    return h * sum;
}

double Integration::simpson_method(std::function<double(double)> f,
                                 double a, double b, int n) {
    if (n <= 0 || a >= b) {
        std::cerr << "Błąd: Niepoprawne parametry całkowania" << std::endl;
        return 0.0;
    }
    
    if (n % 2 != 0) {
        n++; // n musi być parzyste dla metody Simpsona
    }
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x) * (i % 2 == 0 ? 2 : 4);
    }
    
    return h * sum / 3.0;
}

std::vector<Integration::GaussLegendreNode> Integration::get_gauss_legendre_nodes(int n) {
    std::vector<GaussLegendreNode> nodes;
    
    switch (n) {
        case 1:
            nodes.emplace_back(0.0, 2.0);
            break;
            
        case 2:
            nodes.emplace_back(-0.577350269189626, 1.0);
            nodes.emplace_back(0.577350269189626, 1.0);
            break;
            
        case 3:
            nodes.emplace_back(-0.774596669241483, 0.555555555555556);
            nodes.emplace_back(0.0, 0.888888888888889);
            nodes.emplace_back(0.774596669241483, 0.555555555555556);
            break;
            
        case 4:
            nodes.emplace_back(-0.861136311594053, 0.347854845137454);
            nodes.emplace_back(-0.339981043584856, 0.652145154862546);
            nodes.emplace_back(0.339981043584856, 0.652145154862546);
            nodes.emplace_back(0.861136311594053, 0.347854845137454);
            break;
            
        case 5:
            nodes.emplace_back(-0.906179845938664, 0.236926885056189);
            nodes.emplace_back(-0.538469310105683, 0.478628670499366);
            nodes.emplace_back(0.0, 0.568888888888889);
            nodes.emplace_back(0.538469310105683, 0.478628670499366);
            nodes.emplace_back(0.906179845938664, 0.236926885056189);
            break;
            
        default:
            std::cerr << "Błąd: Obsługiwane liczby węzłów G-L: 1-5" << std::endl;
            break;
    }
    
    return nodes;
}

double Integration::gauss_legendre_quadrature(std::function<double(double)> f,
                                            double a, double b, int n) {
    if (a >= b) {
        std::cerr << "Błąd: Niepoprawne granice całkowania" << std::endl;
        return 0.0;
    }
    
    auto nodes = get_gauss_legendre_nodes(n);
    if (nodes.empty()) {
        return 0.0;
    }
    
    double result = 0.0;
    
    for (const auto& node : nodes) {
        // Mapowanie z [-1,1] na [a,b]
        double x = map_to_interval(node.point, a, b);
        double fx = f(x);
        
        if (std::isfinite(fx)) {
            result += node.weight * fx;
        } else {
            std::cerr << "Ostrzeżenie: Nieskończona wartość funkcji w x = " << x << std::endl;
        }
    }
    
    // Mnożenie przez jakobian transformacji
    result *= (b - a) / 2.0;
    
    return result;
}

double Integration::adaptive_gauss_legendre(std::function<double(double)> f,
                                          double a, double b, int n, int num_segments) {
    if (a >= b || num_segments <= 0) {
        return 0.0;
    }
    
    double h = (b - a) / num_segments;
    double sum = 0.0;
    
    for (int i = 0; i < num_segments; ++i) {
        double segment_a = a + i * h;
        double segment_b = a + (i + 1) * h;
        sum += gauss_legendre_quadrature(f, segment_a, segment_b, n);
    }
    
    return sum;
}

Integration::IntegrationResult Integration::adaptive_simpson(
    std::function<double(double)> f,
    double a, double b,
    double tolerance,
    int max_depth) {
    
    IntegrationResult result;
    result.method = "Adaptive Simpson";
    
    // Początkowe obliczenia
    double fa = f(a);
    double fb = f(b);
    double fc = f((a + b) / 2.0);
    double S = simpson_method(f, a, b, 2); // Simpson na całym przedziale
    
    result = adaptive_simpson_recursive(f, a, b, tolerance, S, fa, fb, fc, 
                                      0, max_depth, result.function_calls);
    
    return result;
}

Integration::IntegrationResult Integration::adaptive_simpson_recursive(
    std::function<double(double)> f,
    double a, double b, double tolerance,
    double S, double fa, double fb, double fc,
    int depth, int max_depth, int& function_calls) {
    
    IntegrationResult result;
    
    if (depth >= max_depth) {
        result.value = S;
        result.error_estimate = tolerance;
        result.converged = false;
        return result;
    }
    
    double c = (a + b) / 2.0;
    double fd = f((a + c) / 2.0);
    double fe = f((c + b) / 2.0);
    function_calls += 2;
    
    double S1 = (b - a) / 12.0 * (fa + 4*fd + 2*fc);  // Simpson na [a,c]
    double S2 = (b - a) / 12.0 * (2*fc + 4*fe + fb);  // Simpson na [c,b]
    double S_new = S1 + S2;
    
    double error = std::abs(S_new - S) / 15.0; // Oszacowanie błędu
    
    if (error < tolerance) {
        result.value = S_new + error; // Korekta Richardsona
        result.error_estimate = error;
        result.converged = true;
    } else {
        // Rekursywny podział
        auto left = adaptive_simpson_recursive(f, a, c, tolerance/2, S1, 
                                             fa, fc, fd, depth+1, max_depth, function_calls);
        auto right = adaptive_simpson_recursive(f, c, b, tolerance/2, S2, 
                                              fc, fb, fe, depth+1, max_depth, function_calls);
        
        result.value = left.value + right.value;
        result.error_estimate = left.error_estimate + right.error_estimate;
        result.converged = left.converged && right.converged;
    }
    
    return result;
}

std::vector<Integration::IntegrationResult> Integration::compare_methods(
    std::function<double(double)> f,
    double a, double b, int n,
    double exact_value) {
    
    std::vector<IntegrationResult> results;
    
    // Metoda prostokątów
    IntegrationResult rect_result;
    rect_result.method = "Rectangle";
    rect_result.value = rectangle_method(f, a, b, n);
    rect_result.function_calls = n;
    if (exact_value != 0.0) {
        rect_result.error_estimate = std::abs(rect_result.value - exact_value);
    }
    results.push_back(rect_result);
    
    // Metoda trapezów
    IntegrationResult trap_result;
    trap_result.method = "Trapezoid";
    trap_result.value = trapezoid_method(f, a, b, n);
    trap_result.function_calls = n + 1;
    if (exact_value != 0.0) {
        trap_result.error_estimate = std::abs(trap_result.value - exact_value);
    }
    results.push_back(trap_result);
    
    // Metoda Simpsona
    IntegrationResult simp_result;
    simp_result.method = "Simpson";
    simp_result.value = simpson_method(f, a, b, n);
    simp_result.function_calls = n + 1;
    if (exact_value != 0.0) {
        simp_result.error_estimate = std::abs(simp_result.value - exact_value);
    }
    results.push_back(simp_result);
    
    // Kwadratura Gaussa-Legendre'a (różne liczby węzłów)
    for (int gl_n = 2; gl_n <= 5; ++gl_n) {
        IntegrationResult gl_result;
        gl_result.method = "Gauss-Legendre " + std::to_string(gl_n);
        gl_result.value = gauss_legendre_quadrature(f, a, b, gl_n);
        gl_result.function_calls = gl_n;
        if (exact_value != 0.0) {
            gl_result.error_estimate = std::abs(gl_result.value - exact_value);
        }
        results.push_back(gl_result);
    }
    
    return results;
}

void Integration::test_convergence(std::function<double(double)> f,
                                 double a, double b,
                                 const std::string& method,
                                 double exact_value,
                                 const std::string& filename,
                                 int max_n) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Nie można utworzyć pliku: " << filename << std::endl;
        return;
    }
    
    file << "# Test zbieżności metody: " << method << std::endl;
    file << "# n, integral_value, error, relative_error" << std::endl;
    
    std::cout << "Test zbieżności dla metody: " << method << std::endl;
    
    for (int n = 2; n <= max_n; n *= 2) {
        double result = 0.0;
        
        if (method == "rectangle") {
            result = rectangle_method(f, a, b, n);
        } else if (method == "trapezoid") {
            result = trapezoid_method(f, a, b, n);
        } else if (method == "simpson") {
            result = simpson_method(f, a, b, n);
        } else {
            std::cerr << "Nieznana metoda: " << method << std::endl;
            return;
        }
        
        double error = std::abs(result - exact_value);
        double relative_error = error / std::abs(exact_value);
        
        file << n << ", " << std::scientific << std::setprecision(12)
             << result << ", " << error << ", " << relative_error << std::endl;
        
        std::cout << "n = " << std::setw(6) << n 
                  << ", błąd = " << std::scientific << std::setprecision(3) 
                  << error << std::endl;
    }
    
    file.close();
    std::cout << "Wyniki zapisane do: " << filename << std::endl;
}

double Integration::integrate_oscillatory(std::function<double(double)> f,
                                        double a, double b,
                                        double estimated_period,
                                        int n_per_period) {
    if (estimated_period <= 0 || n_per_period <= 0) {
        return gauss_legendre_quadrature(f, a, b, 5);
    }
    
    // Oblicz liczbę segmentów na podstawie okresu oscylacji
    int num_periods = static_cast<int>(std::ceil((b - a) / estimated_period));
    int num_segments = num_periods * n_per_period;
    
    // Używaj adaptacyjnej kwadratury G-L
    return adaptive_gauss_legendre(f, a, b, 5, num_segments);
}

double Integration::integrate_exponential_growth(std::function<double(double)> f,
                                               double a, double b,
                                               double growth_parameter) {
    if (growth_parameter <= 0) {
        return gauss_legendre_quadrature(f, a, b, 5);
    }
    
    // Dla funkcji exp(ax²) używamy gęstszego podziału w obszarach szybkiego wzrostu
    std::vector<double> breakpoints;
    breakpoints.push_back(a);
    
    // Dodaj więcej punktów w okolicach maksimum
    if (a < 0 && b > 0) {
        // Gęstszy podział wokół zera
        for (double x = -1.0; x <= 1.0; x += 0.2) {
            if (x > a && x < b) {
                breakpoints.push_back(x);
            }
        }
    }
    
    // Dodaj punkty na końcach przedziału
    int num_edge_points = 5;
    double edge_width = 0.1 * (b - a);
    
    for (int i = 1; i < num_edge_points; ++i) {
        double x_left = a + i * edge_width / num_edge_points;
        double x_right = b - i * edge_width / num_edge_points;
        
        if (x_left < b) breakpoints.push_back(x_left);
        if (x_right > a) breakpoints.push_back(x_right);
    }
    
    breakpoints.push_back(b);
    
    // Sortuj punkty podziału
    std::sort(breakpoints.begin(), breakpoints.end());
    breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()), 
                     breakpoints.end());
    
    // Całkuj na każdym podprzedziale
    double total = 0.0;
    for (size_t i = 0; i < breakpoints.size() - 1; ++i) {
        total += gauss_legendre_quadrature(f, breakpoints[i], breakpoints[i+1], 5);
    }
    
    return total;
}

std::vector<std::pair<std::string, double>> Integration::benchmark_methods(
    std::function<double(double)> f,
    double a, double b, int n,
    int iterations) {
    
    std::vector<std::pair<std::string, double>> results;
    
    // Benchmark metody prostokątów
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        volatile double result = rectangle_method(f, a, b, n);
        (void)result;
    }
    auto end = std::chrono::high_resolution_clock::now();
    double rect_time = std::chrono::duration<double, std::milli>(end - start).count();
    results.emplace_back("Rectangle", rect_time);
    
    // Benchmark metody trapezów
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        volatile double result = trapezoid_method(f, a, b, n);
        (void)result;
    }
    end = std::chrono::high_resolution_clock::now();
    double trap_time = std::chrono::duration<double, std::milli>(end - start).count();
    results.emplace_back("Trapezoid", trap_time);
    
    // Benchmark metody Simpsona
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        volatile double result = simpson_method(f, a, b, n);
        (void)result;
    }
    end = std::chrono::high_resolution_clock::now();
    double simp_time = std::chrono::duration<double, std::milli>(end - start).count();
    results.emplace_back("Simpson", simp_time);
    
    // Benchmark kwadratury G-L
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        volatile double result = gauss_legendre_quadrature(f, a, b, 5);
        (void)result;
    }
    end = std::chrono::high_resolution_clock::now();
    double gl_time = std::chrono::duration<double, std::milli>(end - start).count();
    results.emplace_back("Gauss-Legendre", gl_time);
    
    return results;
}

bool Integration::run_validation_tests() {
    const double tolerance = 1e-10;
    bool all_passed = true;
    
    std::cout << "Uruchamianie testów walidacyjnych metod całkowania..." << std::endl;
    
    // Test 1: Wielomian stopnia 2 - dokładny dla metody Simpsona
    auto f1 = [](double x) { return x*x; };
    double exact1 = 8.0/3.0; // ∫₀² x² dx = 8/3
    double result1 = simpson_method(f1, 0.0, 2.0, 2);
    
    if (std::abs(result1 - exact1) < tolerance) {
        std::cout << "✓ Test wielomianu x² - POWODZENIE" << std::endl;
    } else {
        std::cout << "✗ Test wielomianu x² - BŁĄD (oczekiwane: " << exact1 
                  << ", otrzymane: " << result1 << ")" << std::endl;
        all_passed = false;
    }
    
    // Test 2: Funkcja stała - dokładna dla wszystkich metod
    auto f2 = [](double x) { (void)x; return 5.0; };
    double exact2 = 10.0; // ∫₁³ 5 dx = 10
    double result2_rect = rectangle_method(f2, 1.0, 3.0, 10);
    double result2_trap = trapezoid_method(f2, 1.0, 3.0, 10);
    double result2_simp = simpson_method(f2, 1.0, 3.0, 10);
    
    if (std::abs(result2_rect - exact2) < tolerance && 
        std::abs(result2_trap - exact2) < tolerance &&
        std::abs(result2_simp - exact2) < tolerance) {
        std::cout << "✓ Test funkcji stałej - POWODZENIE" << std::endl;
    } else {
        std::cout << "✗ Test funkcji stałej - BŁĄD" << std::endl;
        all_passed = false;
    }
    
    // Test 3: Kwadratura G-L dla wielomianu stopnia 3 (4 węzły powinny być dokładne)
    auto f3 = [](double x) { return x*x*x; };
    double exact3 = 0.0; // ∫₋₁¹ x³ dx = 0 (funkcja nieparzysta)
    double result3 = gauss_legendre_quadrature(f3, -1.0, 1.0, 4);
    
    if (std::abs(result3 - exact3) < tolerance) {
        std::cout << "✓ Test kwadratury G-L dla x³ - POWODZENIE" << std::endl;
    } else {
        std::cout << "✗ Test kwadratury G-L dla x³ - BŁĄD (otrzymane: " 
                  << result3 << ")" << std::endl;
        all_passed = false;
    }
    
    // Test 4: Funkcja exp(x)
    auto f4 = [](double x) { return std::exp(x); };
    double exact4 = std::exp(1.0) - 1.0; // ∫₀¹ eˣ dx = e - 1
    double result4 = adaptive_simpson(f4, 0.0, 1.0, 1e-8, 10).value;
    
    if (std::abs(result4 - exact4) < 1e-6) {
        std::cout << "✓ Test adaptacyjny Simpson dla exp(x) - POWODZENIE" << std::endl;
    } else {
        std::cout << "✗ Test adaptacyjny Simpson dla exp(x) - BŁĄD" << std::endl;
        all_passed = false;
    }
    
    return all_passed;
}

} // namespace agh_numerical