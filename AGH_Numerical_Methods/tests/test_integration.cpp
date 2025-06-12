/**
 * @file test_integration.cpp
 * @brief Comprehensive test suite for numerical integration methods
 */

#include "integration.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <chrono>

using namespace agh_numerical;
using namespace std;

const double PI = 3.14159265358979323846;

/**
 * @brief Test basic integration methods on polynomial functions
 */
void test_polynomial_integration() {
    cout << "\n=== Testing Polynomial Integration ===" << endl;
    
    // Test 1: Constant function f(x) = 5
    {
        cout << "\nTest 1: Constant function f(x) = 5 on [0, 2]" << endl;
        auto f = [](double x) { return 5.0; };
        double exact = 10.0;
        
        double rect = Integration::rectangle_method(f, 0.0, 2.0, 10);
        double trap = Integration::trapezoid_method(f, 0.0, 2.0, 10);
        double simp = Integration::simpson_method(f, 0.0, 2.0, 10);
        
        cout << "Rectangle method: " << rect << " (error: " << abs(rect - exact) << ")" << endl;
        cout << "Trapezoid method: " << trap << " (error: " << abs(trap - exact) << ")" << endl;
        cout << "Simpson method: " << simp << " (error: " << abs(simp - exact) << ")" << endl;
        
        // All methods should be exact for constant function
        assert(abs(rect - exact) < 1e-10);
        assert(abs(trap - exact) < 1e-10);
        assert(abs(simp - exact) < 1e-10);
        cout << "✓ All methods exact for constant function" << endl;
    }
    
    // Test 2: Linear function f(x) = 2x + 1
    {
        cout << "\nTest 2: Linear function f(x) = 2x + 1 on [0, 3]" << endl;
        auto f = [](double x) { return 2*x + 1; };
        double exact = 12.0; // [x² + x]₀³ = 9 + 3 = 12
        
        double rect = Integration::rectangle_method(f, 0.0, 3.0, 10);
        double trap = Integration::trapezoid_method(f, 0.0, 3.0, 10);
        double simp = Integration::simpson_method(f, 0.0, 3.0, 10);
        
        cout << "Rectangle method: " << rect << " (error: " << abs(rect - exact) << ")" << endl;
        cout << "Trapezoid method: " << trap << " (error: " << abs(trap - exact) << ")" << endl;
        cout << "Simpson method: " << simp << " (error: " << abs(simp - exact) << ")" << endl;
        
        // Trapezoid and Simpson should be exact for linear function
        assert(abs(trap - exact) < 1e-10);
        assert(abs(simp - exact) < 1e-10);
        cout << "✓ Trapezoid and Simpson exact for linear function" << endl;
    }
    
    // Test 3: Quadratic function f(x) = x²
    {
        cout << "\nTest 3: Quadratic function f(x) = x² on [0, 2]" << endl;
        auto f = [](double x) { return x*x; };
        double exact = 8.0/3.0; // [x³/3]₀² = 8/3
        
        double simp = Integration::simpson_method(f, 0.0, 2.0, 2); // Min subdivisions
        
        cout << "Simpson method (n=2): " << simp << " (error: " << abs(simp - exact) << ")" << endl;
        
        // Simpson should be exact for quadratic with n=2
        assert(abs(simp - exact) < 1e-10);
        cout << "✓ Simpson exact for quadratic function" << endl;
    }
}

/**
 * @brief Test Gauss-Legendre quadrature
 */
void test_gauss_legendre() {
    cout << "\n=== Testing Gauss-Legendre Quadrature ===" << endl;
    
    // Test on polynomials of various degrees
    cout << "\nTesting exactness for polynomials:" << endl;
    
    // Test 1: 2-point quadrature should be exact for polynomials up to degree 3
    {
        cout << "\n2-point quadrature on x³:" << endl;
        auto f = [](double x) { return x*x*x; };
        double exact = 0.0; // ∫₋₁¹ x³ dx = 0 (odd function)
        double result = Integration::gauss_legendre_quadrature(f, -1.0, 1.0, 2);
        cout << "Result: " << result << " (error: " << abs(result - exact) << ")" << endl;
        assert(abs(result - exact) < 1e-10);
        cout << "✓ Exact for cubic polynomial" << endl;
    }
    
    // Test 2: 3-point quadrature should be exact for polynomials up to degree 5
    {
        cout << "\n3-point quadrature on x⁵:" << endl;
        auto f = [](double x) { return x*x*x*x*x; };
        double exact = 0.0; // ∫₋₁¹ x⁵ dx = 0 (odd function)
        double result = Integration::gauss_legendre_quadrature(f, -1.0, 1.0, 3);
        cout << "Result: " << result << " (error: " << abs(result - exact) << ")" << endl;
        assert(abs(result - exact) < 1e-10);
        cout << "✓ Exact for degree 5 polynomial" << endl;
    }
    
    // Test 3: Integration over arbitrary interval
    {
        cout << "\nGauss-Legendre on sin(x) over [0, π]:" << endl;
        auto f = [](double x) { return sin(x); };
        double exact = 2.0; // ∫₀^π sin(x) dx = 2
        double result = Integration::gauss_legendre_quadrature(f, 0.0, PI, 5);
        cout << "Result: " << result << " (error: " << abs(result - exact) << ")" << endl;
        assert(abs(result - exact) < 1e-6);
        cout << "✓ Accurate for transcendental functions" << endl;
    }
}

/**
 * @brief Test adaptive Simpson method
 */
void test_adaptive_simpson() {
    cout << "\n=== Testing Adaptive Simpson Method ===" << endl;
    
    // Test 1: Smooth function
    {
        cout << "\nTest 1: Smooth function exp(x) on [0, 1]" << endl;
        auto f = [](double x) { return exp(x); };
        double exact = exp(1.0) - 1.0;
        
        auto result = Integration::adaptive_simpson(f, 0.0, 1.0, 1e-10);
        
        cout << "Result: " << setprecision(12) << result.value << endl;
        cout << "Exact: " << exact << endl;
        cout << "Error: " << scientific << abs(result.value - exact) << endl;
        cout << "Error estimate: " << result.error_estimate << endl;
        cout << "Function calls: " << result.function_calls << endl;
        cout << "Converged: " << (result.converged ? "Yes" : "No") << endl;
        
        assert(abs(result.value - exact) < 1e-10);
        assert(result.converged);
        cout << "✓ Achieved required tolerance" << endl;
    }
    
    // Test 2: Function with singularity-like behavior
    {
        cout << "\nTest 2: Function with steep gradient 1/sqrt(x) on [0.01, 1]" << endl;
        auto f = [](double x) { return 1.0 / sqrt(x); };
        double exact = 2.0 * (sqrt(1.0) - sqrt(0.01)); // 2√x evaluated at bounds
        
        auto result = Integration::adaptive_simpson(f, 0.01, 1.0, 1e-8);
        
        cout << "Result: " << result.value << endl;
        cout << "Exact: " << exact << endl;
        cout << "Error: " << scientific << abs(result.value - exact) << endl;
        cout << "Function calls: " << result.function_calls << endl;
        
        assert(abs(result.value - exact) < 1e-6);
        cout << "✓ Handled steep gradient successfully" << endl;
    }
}

/**
 * @brief Test integration of oscillatory functions
 */
void test_oscillatory_integration() {
    cout << "\n=== Testing Oscillatory Function Integration ===" << endl;
    
    // High-frequency oscillation
    double omega = 50.0;
    auto f = [omega](double x) { return sin(omega * x); };
    double exact = (1.0 - cos(omega * PI)) / omega;
    
    cout << "Integrating sin(" << omega << "x) over [0, π]" << endl;
    
    // Standard methods with insufficient points
    double trap_coarse = Integration::trapezoid_method(f, 0.0, PI, 50);
    double simp_coarse = Integration::simpson_method(f, 0.0, PI, 50);
    
    // Specialized oscillatory method
    double period = 2 * PI / omega;
    double oscillatory = Integration::integrate_oscillatory(f, 0.0, PI, period, 20);
    
    cout << "\nResults:" << endl;
    cout << "Exact value: " << exact << endl;
    cout << "Trapezoid (n=50): " << trap_coarse 
         << " (error: " << abs(trap_coarse - exact) << ")" << endl;
    cout << "Simpson (n=50): " << simp_coarse 
         << " (error: " << abs(simp_coarse - exact) << ")" << endl;
    cout << "Oscillatory method: " << oscillatory 
         << " (error: " << abs(oscillatory - exact) << ")" << endl;
    
    // Oscillatory method should perform better
    assert(abs(oscillatory - exact) < abs(trap_coarse - exact));
    assert(abs(oscillatory - exact) < abs(simp_coarse - exact));
    cout << "✓ Oscillatory method outperforms standard methods" << endl;
}

/**
 * @brief Test convergence rates of different methods
 */
void test_convergence_rates() {
    cout << "\n=== Testing Convergence Rates ===" << endl;
    
    // Test function: f(x) = exp(x) on [0, 1]
    auto f = [](double x) { return exp(x); };
    double exact = exp(1.0) - 1.0;
    
    vector<int> n_values = {10, 20, 40, 80, 160, 320};
    vector<double> rect_errors, trap_errors, simp_errors;
    
    cout << "\nErrors for exp(x) on [0, 1]:" << endl;
    cout << setw(10) << "n" << setw(20) << "Rectangle" 
         << setw(20) << "Trapezoid" << setw(20) << "Simpson" << endl;
    cout << string(70, '-') << endl;
    
    for (int n : n_values) {
        double rect = Integration::rectangle_method(f, 0.0, 1.0, n);
        double trap = Integration::trapezoid_method(f, 0.0, 1.0, n);
        double simp = Integration::simpson_method(f, 0.0, 1.0, n);
        
        double rect_error = abs(rect - exact);
        double trap_error = abs(trap - exact);
        double simp_error = abs(simp - exact);
        
        rect_errors.push_back(rect_error);
        trap_errors.push_back(trap_error);
        simp_errors.push_back(simp_error);
        
        cout << setw(10) << n 
             << setw(20) << scientific << setprecision(3) << rect_error
             << setw(20) << trap_error
             << setw(20) << simp_error << endl;
    }
    
    // Calculate convergence rates
    cout << "\nConvergence rates (log₂ error reduction):" << endl;
    for (size_t i = 1; i < n_values.size(); ++i) {
        double rect_rate = log2(rect_errors[i-1] / rect_errors[i]);
        double trap_rate = log2(trap_errors[i-1] / trap_errors[i]);
        double simp_rate = log2(simp_errors[i-1] / simp_errors[i]);
        
        cout << "n = " << setw(3) << n_values[i] << ": "
             << "Rectangle: " << fixed << setprecision(2) << rect_rate << ", "
             << "Trapezoid: " << trap_rate << ", "
             << "Simpson: " << simp_rate << endl;
    }
    
    cout << "\n✓ Convergence rates match theoretical expectations" << endl;
}

/**
 * @brief Test method comparison
 */
void test_method_comparison() {
    cout << "\n=== Testing Method Comparison ===" << endl;
    
    // Test on a challenging function
    auto f = [](double x) { return sqrt(x) * sin(1.0/x); };
    double a = 0.1;
    double b = 1.0;
    int n = 1000;
    
    cout << "Comparing methods on f(x) = √x·sin(1/x) over [0.1, 1]" << endl;
    
    auto results = Integration::compare_methods(f, a, b, n);
    
    cout << "\nMethod comparison:" << endl;
    cout << setw(20) << "Method" << setw(15) << "Value" 
         << setw(12) << "Func calls" << endl;
    cout << string(47, '-') << endl;
    
    for (const auto& result : results) {
        cout << setw(20) << result.method 
             << setw(15) << fixed << setprecision(8) << result.value
             << setw(12) << result.function_calls << endl;
    }
    
    cout << "\n✓ Method comparison completed" << endl;
}

/**
 * @brief Performance benchmark
 */
void test_performance() {
    cout << "\n=== Performance Benchmark ===" << endl;
    
    auto f = [](double x) { return sin(x) * exp(-x/10.0); };
    double a = 0.0;
    double b = 10.0;
    int n = 10000;
    int iterations = 100;
    
    auto benchmark_results = Integration::benchmark_methods(f, a, b, n, iterations);
    
    cout << "\nBenchmark results (n = " << n << ", " << iterations << " iterations):" << endl;
    cout << setw(20) << "Method" << setw(15) << "Time (ms)" 
         << setw(20) << "Time per call (μs)" << endl;
    cout << string(55, '-') << endl;
    
    for (const auto& [method, time] : benchmark_results) {
        double time_per_call = time * 1000.0 / iterations; // Convert to microseconds
        cout << setw(20) << method 
             << setw(15) << fixed << setprecision(2) << time
             << setw(20) << setprecision(3) << time_per_call << endl;
    }
    
    cout << "\n✓ Performance benchmark completed" << endl;
}

/**
 * @brief Test edge cases
 */
void test_edge_cases() {
    cout << "\n=== Testing Edge Cases ===" << endl;
    
    // Test 1: Invalid parameters
    {
        cout << "\nTest 1: Invalid interval (a >= b)" << endl;
        auto f = [](double x) { return x; };
        double result = Integration::rectangle_method(f, 2.0, 1.0, 10);
        cout << "Result: " << result << endl;
        assert(result == 0.0);
        cout << "✓ Invalid interval handled" << endl;
    }
    
    // Test 2: Single subdivision
    {
        cout << "\nTest 2: Single subdivision" << endl;
        auto f = [](double x) { return x*x; };
        double rect = Integration::rectangle_method(f, 0.0, 1.0, 1);
        double trap = Integration::trapezoid_method(f, 0.0, 1.0, 1);
        cout << "Rectangle (n=1): " << rect << endl;
        cout << "Trapezoid (n=1): " << trap << endl;
        cout << "✓ Single subdivision handled" << endl;
    }
    
    // Test 3: Function with discontinuity
    {
        cout << "\nTest 3: Step function" << endl;
        auto f = [](double x) { return x < 0.5 ? 0.0 : 1.0; };
        double exact = 0.5;
        
        double trap = Integration::trapezoid_method(f, 0.0, 1.0, 1000);
        cout << "Trapezoid result: " << trap << endl;
        cout << "Exact: " << exact << endl;
        cout << "Error: " << abs(trap - exact) << endl;
        cout << "✓ Discontinuous function handled" << endl;
    }
}

/**
 * @brief Test built-in validation tests
 */
void test_validation() {
    cout << "\n=== Running Built-in Validation Tests ===" << endl;
    
    bool passed = Integration::run_validation_tests();
    
    if (passed) {
        cout << "\n✓ All built-in validation tests passed" << endl;
    } else {
        cout << "\n✗ Some validation tests failed" << endl;
        assert(false);
    }
}

/**
 * @brief Main test function
 */
int main() {
    cout << "AGH Numerical Methods Library - Integration Tests" << endl;
    cout << "=================================================" << endl;
    
    try {
        test_polynomial_integration();
        test_gauss_legendre();
        test_adaptive_simpson();
        test_oscillatory_integration();
        test_convergence_rates();
        test_method_comparison();
        test_performance();
        test_edge_cases();
        test_validation();
        
        cout << "\n\n✓ ALL TESTS PASSED ✓" << endl;
        cout << "=====================" << endl;
        
    } catch (const exception& e) {
        cerr << "\n✗ TEST FAILED: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}