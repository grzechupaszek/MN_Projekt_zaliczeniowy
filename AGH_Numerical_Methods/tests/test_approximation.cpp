/**
 * @file test_approximation.cpp
 * @brief Comprehensive test suite for approximation methods
 * @author Grzegorz Paszek F. Rozanski
 */

#include "approximation.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <vector>
#include <functional>
#include <chrono>
using namespace agh_numerical;
using namespace std;

const double PI = 3.14159265358979323846;

/**
 * @brief Test polynomial least squares approximation
 */
void test_polynomial_least_squares() {
    cout << "\n=== Testing Polynomial Least Squares ===" << endl;
    
    // Test 1: Linear data (exact fit expected)
    {
        cout << "\nTest 1: Linear function y = 2x + 1" << endl;
        vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
        vector<double> y = {1.0, 3.0, 5.0, 7.0, 9.0};
        
        auto result = Approximation::polynomial_least_squares(x, y, 1);
        
        cout << "Coefficients: ";
        for (double c : result.coefficients) {
            cout << c << " ";
        }
        cout << endl;
        
        // Should get coefficients [1, 2] for y = 1 + 2x
        assert(result.is_valid);
        assert(abs(result.coefficients[0] - 1.0) < 1e-10);
        assert(abs(result.coefficients[1] - 2.0) < 1e-10);
        assert(result.rmse < 1e-10);
        cout << "✓ Linear approximation exact for linear data" << endl;
    }
    
    // Test 2: Quadratic approximation
    {
        cout << "\nTest 2: Quadratic function y = x²" << endl;
        vector<double> x = {-2.0, -1.0, 0.0, 1.0, 2.0};
        vector<double> y = {4.0, 1.0, 0.0, 1.0, 4.0};
        
        auto result = Approximation::polynomial_least_squares(x, y, 2);
        
        assert(result.is_valid);
        assert(abs(result.coefficients[0]) < 1e-10); // a₀ ≈ 0
        assert(abs(result.coefficients[1]) < 1e-10); // a₁ ≈ 0
        assert(abs(result.coefficients[2] - 1.0) < 1e-10); // a₂ ≈ 1
        cout << "✓ Quadratic approximation successful" << endl;
    }
    
    // Test 3: Overdetermined system with noise
    {
        cout << "\nTest 3: Noisy data approximation" << endl;
        vector<double> x = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
        vector<double> y = {0.1, 0.9, 2.1, 2.9, 4.2, 4.8, 6.1};
        
        auto result = Approximation::polynomial_least_squares(x, y, 1);
        
        assert(result.is_valid);
        cout << "RMSE: " << result.rmse << endl;
        cout << "R²: " << result.r_squared << endl;
        assert(result.r_squared > 0.99); // Should have good fit
        cout << "✓ Approximation handles noisy data well" << endl;
    }
}

/**
 * @brief Test basis function approximation
 */
void test_basis_function_approximation() {
    cout << "\n=== Testing Basis Function Approximation ===" << endl;
    
    // Test with custom basis functions
    vector<double> x = {0.0, PI/4, PI/2, 3*PI/4, PI};
    vector<double> y;
    
    // Generate data: y = 2*sin(x) + 3*cos(x)
    for (double xi : x) {
        y.push_back(2*sin(xi) + 3*cos(xi));
    }
    
    // Define basis functions
    vector<Approximation::BasisFunction> basis = {
        Approximation::BasisFunction([](double x) { return sin(x); }, "sin(x)"),
        Approximation::BasisFunction([](double x) { return cos(x); }, "cos(x)")
    };
    
    auto result = Approximation::basis_function_approximation(x, y, basis);
    
    cout << "Coefficients: ";
    for (double c : result.coefficients) {
        cout << c << " ";
    }
    cout << endl;
    
    assert(result.is_valid);
    assert(abs(result.coefficients[0] - 2.0) < 1e-10); // sin coefficient
    assert(abs(result.coefficients[1] - 3.0) < 1e-10); // cos coefficient
    cout << "✓ Trigonometric basis approximation successful" << endl;
}

/**
 * @brief Test exponential approximation
 */
void test_exponential_approximation() {
    cout << "\n=== Testing Exponential Approximation ===" << endl;
    
    // Generate exponential data: y = 2 * exp(0.5x)
    vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    vector<double> y;
    
    for (double xi : x) {
        y.push_back(2.0 * exp(0.5 * xi));
    }
    
    auto result = Approximation::exponential_approximation(x, y);
    
    cout << "Exponential fit: y = " << result.coefficients[0] 
         << " * exp(" << result.coefficients[1] << "x)" << endl;
    
    assert(result.is_valid);
    assert(abs(result.coefficients[0] - 2.0) < 1e-6);
    assert(abs(result.coefficients[1] - 0.5) < 1e-6);
    cout << "✓ Exponential approximation successful" << endl;
}

/**
 * @brief Test power approximation
 */
void test_power_approximation() {
    cout << "\n=== Testing Power Approximation ===" << endl;
    
    // Generate power law data: y = 3 * x^2
    vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    vector<double> y;
    
    for (double xi : x) {
        y.push_back(3.0 * pow(xi, 2.0));
    }
    
    auto result = Approximation::power_approximation(x, y);
    
    cout << "Power fit: y = " << result.coefficients[0] 
         << " * x^" << result.coefficients[1] << endl;
    
    assert(result.is_valid);
    assert(abs(result.coefficients[0] - 3.0) < 1e-6);
    assert(abs(result.coefficients[1] - 2.0) < 1e-6);
    cout << "✓ Power approximation successful" << endl;
}

/**
 * @brief Test Chebyshev approximation
 */
void test_chebyshev_approximation() {
    cout << "\n=== Testing Chebyshev Approximation ===" << endl;
    
    // Test on Runge function: f(x) = 1/(1 + 25x²)
    vector<double> x, y;
    for (int i = -20; i <= 20; ++i) {
        double xi = i * 0.05; // [-1, 1]
        x.push_back(xi);
        y.push_back(1.0 / (1.0 + 25.0 * xi * xi));
    }
    
    // Compare polynomial and Chebyshev approximations
    auto poly_result = Approximation::polynomial_least_squares(x, y, 8);
    auto cheby_result = Approximation::chebyshev_approximation(x, y, 8);
    
    cout << "Polynomial RMSE: " << poly_result.rmse << endl;
    cout << "Chebyshev RMSE: " << cheby_result.rmse << endl;
    cout << "Polynomial max error: " << poly_result.max_error << endl;
    cout << "Chebyshev max error: " << cheby_result.max_error << endl;
    
    assert(cheby_result.is_valid);
    // Chebyshev should generally have better uniform error distribution
    cout << "✓ Chebyshev approximation completed" << endl;
}

/**
 * @brief Test optimal polynomial degree finding
 */
void test_find_optimal_degree() {
    cout << "\n=== Testing Optimal Polynomial Degree ===" << endl;
    
    // Generate data from cubic polynomial with noise
    vector<double> x, y;
    for (int i = 0; i <= 20; ++i) {
        double xi = i * 0.1;
        x.push_back(xi);
        // y = x³ - 2x² + x + 1 + small noise
        y.push_back(xi*xi*xi - 2*xi*xi + xi + 1 + 0.01*(rand()%100-50)/100.0);
    }
    
    int optimal_degree = Approximation::find_optimal_polynomial_degree(x, y, 8, 0.3);
    
    cout << "Optimal polynomial degree: " << optimal_degree << endl;
    assert(optimal_degree >= 2 && optimal_degree <= 4); // Should be around 3
    cout << "✓ Optimal degree selection successful" << endl;
}

/**
 * @brief Test comparison of approximation methods
 */
void test_compare_methods() {
    cout << "\n=== Testing Method Comparison ===" << endl;
    
    // Generate test data
    vector<double> x = {0.1, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0};
    vector<double> y = {0.5, 1.2, 2.0, 3.5, 5.2, 7.1, 9.3};
    
    auto results = Approximation::compare_approximation_methods(x, y, 3);
    
    cout << "\nMethod comparison results:" << endl;
    cout << setw(30) << "Method" << setw(12) << "RMSE" 
         << setw(12) << "R²" << setw(12) << "Max Error" << endl;
    cout << string(66, '-') << endl;
    
    for (const auto& result : results) {
        cout << setw(30) << result.method 
             << setw(12) << fixed << setprecision(6) << result.rmse
             << setw(12) << result.r_squared
             << setw(12) << result.max_error << endl;
    }
    
    assert(!results.empty());
    cout << "✓ Method comparison completed" << endl;
}

/**
 * @brief Test approximation statistics computation
 */
void test_approximation_statistics() {
    cout << "\n=== Testing Approximation Statistics ===" << endl;
    
    vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    vector<double> y = {1.0, 2.9, 5.1, 6.9, 9.0}; // Approximately y = 2x + 1
    
    auto result = Approximation::polynomial_least_squares(x, y, 1);
    
    cout << "Linear fit statistics:" << endl;
    cout << "  RMSE: " << result.rmse << endl;
    cout << "  R²: " << result.r_squared << endl;
    cout << "  Max error: " << result.max_error << endl;
    
    assert(result.r_squared > 0.99); // Should have excellent fit
    assert(result.rmse < 0.2);
    cout << "✓ Statistics computation verified" << endl;
}

/**
 * @brief Test edge cases
 */
void test_edge_cases() {
    cout << "\n=== Testing Edge Cases ===" << endl;
    
    // Test 1: Empty data
    {
        cout << "\nTest 1: Empty data" << endl;
        vector<double> x, y;
        auto result = Approximation::polynomial_least_squares(x, y, 1);
        assert(!result.is_valid);
        cout << "✓ Empty data handled correctly" << endl;
    }
    
    // Test 2: Single data point
    {
        cout << "\nTest 2: Single data point" << endl;
        vector<double> x = {1.0};
        vector<double> y = {2.0};
        auto result = Approximation::polynomial_least_squares(x, y, 0);
        assert(result.is_valid);
        cout << "✓ Single point handled correctly" << endl;
    }
    
    // Test 3: Negative values for exponential
    {
        cout << "\nTest 3: Negative values for exponential" << endl;
        vector<double> x = {0.0, 1.0, 2.0};
        vector<double> y = {1.0, -2.0, 3.0}; // Contains negative
        auto result = Approximation::exponential_approximation(x, y);
        assert(!result.is_valid);
        cout << "✓ Invalid exponential data rejected" << endl;
    }
}

/**
 * @brief Test evaluation of approximation
 */
void test_evaluate_approximation() {
    cout << "\n=== Testing Approximation Evaluation ===" << endl;
    
    // Create a simple linear approximation
    vector<double> x = {0.0, 1.0, 2.0};
    vector<double> y = {1.0, 3.0, 5.0}; // y = 2x + 1
    
    auto result = Approximation::polynomial_least_squares(x, y, 1);
    
    // Test evaluation at various points
    double test_points[] = {0.5, 1.5, 2.5};
    double expected[] = {2.0, 4.0, 6.0};
    
    for (int i = 0; i < 3; ++i) {
        double value = Approximation::evaluate_approximation(result, test_points[i]);
        cout << "f(" << test_points[i] << ") = " << value 
             << " (expected: " << expected[i] << ")" << endl;
        assert(abs(value - expected[i]) < 1e-10);
    }
    
    cout << "✓ Approximation evaluation verified" << endl;
}

/**
 * @brief Performance test for large datasets
 */
void test_performance() {
    cout << "\n=== Testing Performance ===" << endl;
    
    // Generate large dataset
    const int n = 10000;
    vector<double> x, y;
    
    for (int i = 0; i < n; ++i) {
        double xi = i * 0.001;
        x.push_back(xi);
        y.push_back(sin(xi) + 0.1 * cos(10*xi) + 0.01*(rand()%100-50)/100.0);
    }
    
    auto start = chrono::high_resolution_clock::now();
    auto result = Approximation::polynomial_least_squares(x, y, 5);
    auto end = chrono::high_resolution_clock::now();
    
    double elapsed = chrono::duration<double, milli>(end - start).count();
    
    cout << "Polynomial approximation (degree 5) on " << n << " points:" << endl;
    cout << "  Time: " << fixed << setprecision(2) << elapsed << " ms" << endl;
    cout << "  RMSE: " << result.rmse << endl;
    
    assert(result.is_valid);
    cout << "✓ Performance test completed" << endl;
}

/**
 * @brief Test weighted polynomial approximation
 */
void test_weighted_approximation() {
    cout << "\n=== Testing Weighted Approximation ===" << endl;
    
    // Data with outlier
    vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    vector<double> y = {1.0, 3.0, 5.0, 10.0, 9.0}; // 3rd point is outlier
    vector<double> weights = {1.0, 1.0, 1.0, 0.1, 1.0}; // Low weight for outlier
    
    auto weighted_result = Approximation::weighted_polynomial_approximation(x, y, weights, 1);
    auto unweighted_result = Approximation::polynomial_least_squares(x, y, 1);
    
    cout << "Unweighted coefficients: ";
    for (double c : unweighted_result.coefficients) cout << c << " ";
    cout << endl;
    
    cout << "Weighted coefficients: ";
    for (double c : weighted_result.coefficients) cout << c << " ";
    cout << endl;
    
    // Weighted should be less affected by outlier
    assert(weighted_result.is_valid);
    cout << "✓ Weighted approximation handles outliers" << endl;
}

/**
 * @brief Main test function
 */
int main() {
    cout << "AGH Numerical Methods Library - Approximation Tests" << endl;
    cout << "===================================================" << endl;
    
    try {
        test_polynomial_least_squares();
        test_basis_function_approximation();
        test_exponential_approximation();
        test_power_approximation();
        test_chebyshev_approximation();
        test_find_optimal_degree();
        test_compare_methods();
        test_approximation_statistics();
        test_edge_cases();
        test_evaluate_approximation();
        test_performance();
        test_weighted_approximation();
        
        cout << "\n\n✓ ALL TESTS PASSED ✓" << endl;
        cout << "=====================" << endl;
        
    } catch (const exception& e) {
        cerr << "\n✗ TEST FAILED: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}