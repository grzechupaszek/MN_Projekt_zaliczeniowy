/**
 * @file approximation.cpp
 * @brief Implementacja metod aproksymacji funkcji
 */

#include "approximation.h"
#include "linear_systems.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace agh_numerical {

const double Approximation::EPSILON = 1e-14;

Approximation::ApproximationResult Approximation::polynomial_least_squares(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    int degree) {
    
    ApproximationResult result;
    result.method = "Polynomial Least Squares";
    result.degree = degree;
    
    if (x_data.size() != y_data.size() || x_data.empty()) {
        std::cerr << "Błąd: Niezgodne rozmiary danych lub puste dane" << std::endl;
        return result;
    }
    
    if (degree < 0 || degree >= static_cast<int>(x_data.size())) {
        std::cerr << "Błąd: Niepoprawny stopień wielomianu" << std::endl;
        return result;
    }
    
    int n = x_data.size();
    int m = degree + 1; // Liczba współczynników
    
    // Budowa macierzy A układu równań normalnych A^T * A * c = A^T * y
    std::vector<std::vector<double>> ATA(m, std::vector<double>(m, 0.0));
    std::vector<double> ATy(m, 0.0);
    
    // Wypełnienie macierzy A^T * A i wektora A^T * y
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < n; ++k) {
                ATA[i][j] += std::pow(x_data[k], i + j);
            }
        }
        
        for (int k = 0; k < n; ++k) {
            ATy[i] += y_data[k] * std::pow(x_data[k], i);
        }
    }
    
    // Rozwiązanie układu równań normalnych
    auto solution = LinearSystems::gauss_elimination(ATA, ATy);
    
    if (!solution.is_valid) {
        std::cerr << "Błąd: Nie można rozwiązać układu równań normalnych" << std::endl;
        return result;
    }
    
    result.coefficients = solution.x;
    result.is_valid = true;
    
    // Obliczenie statystyk jakości aproksymacji
    result = compute_approximation_statistics(x_data, y_data, result);
    
    return result;
}

Approximation::ApproximationResult Approximation::basis_function_approximation(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<BasisFunction>& basis_functions) {
    
    ApproximationResult result;
    result.method = "Basis Function Approximation";
    
    if (x_data.size() != y_data.size() || x_data.empty() || basis_functions.empty()) {
        return result;
    }
    
    int n = x_data.size();
    int m = basis_functions.size();
    
    // Budowa macierzy układu równań normalnych
    std::vector<std::vector<double>> ATA(m, std::vector<double>(m, 0.0));
    std::vector<double> ATy(m, 0.0);
    
    // Obliczenie wartości funkcji bazowych dla wszystkich punktów
    std::vector<std::vector<double>> phi(n, std::vector<double>(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            phi[i][j] = basis_functions[j].function(x_data[i]);
        }
    }
    
    // Wypełnienie macierzy A^T * A i wektora A^T * y
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < n; ++k) {
                ATA[i][j] += phi[k][i] * phi[k][j];
            }
        }
        
        for (int k = 0; k < n; ++k) {
            ATy[i] += y_data[k] * phi[k][i];
        }
    }
    
    // Rozwiązanie układu
    auto solution = LinearSystems::gauss_elimination(ATA, ATy);
    
    if (solution.is_valid) {
        result.coefficients = solution.x;
        result.is_valid = true;
        result = compute_approximation_statistics(x_data, y_data, result);
    }
    
    return result;
}

Approximation::ApproximationResult Approximation::exponential_approximation(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data) {
    
    ApproximationResult result;
    result.method = "Exponential Approximation";
    
    if (x_data.size() != y_data.size() || x_data.empty()) {
        return result;
    }
    
    // Sprawdzenie czy wszystkie y > 0 (wymagane dla ln(y))
    for (double y : y_data) {
        if (y <= 0) {
            std::cerr << "Błąd: Aproksymacja eksponencjalna wymaga y > 0" << std::endl;
            return result;
        }
    }
    
    // Linearyzacja: ln(y) = ln(a) + bx
    std::vector<double> ln_y(y_data.size());
    for (size_t i = 0; i < y_data.size(); ++i) {
        ln_y[i] = std::log(y_data[i]);
    }
    
    // Regresja liniowa dla ln(y) = c0 + c1*x
    auto linear_result = polynomial_least_squares(x_data, ln_y, 1);
    
    if (linear_result.is_valid) {
        // Przekształcenie z powrotem: a = exp(c0), b = c1
        result.coefficients.resize(2);
        result.coefficients[0] = std::exp(linear_result.coefficients[0]); // a
        result.coefficients[1] = linear_result.coefficients[1];           // b
        result.is_valid = true;
        result = compute_approximation_statistics(x_data, y_data, result);
    }
    
    return result;
}

Approximation::ApproximationResult Approximation::power_approximation(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data) {
    
    ApproximationResult result;
    result.method = "Power Approximation";
    
    if (x_data.size() != y_data.size() || x_data.empty()) {
        return result;
    }
    
    // Sprawdzenie czy wszystkie x > 0 i y > 0
    for (size_t i = 0; i < x_data.size(); ++i) {
        if (x_data[i] <= 0 || y_data[i] <= 0) {
            std::cerr << "Błąd: Aproksymacja potęgowa wymaga x > 0 i y > 0" << std::endl;
            return result;
        }
    }
    
    // Linearyzacja: ln(y) = ln(a) + b*ln(x)
    std::vector<double> ln_x(x_data.size());
    std::vector<double> ln_y(y_data.size());
    
    for (size_t i = 0; i < x_data.size(); ++i) {
        ln_x[i] = std::log(x_data[i]);
        ln_y[i] = std::log(y_data[i]);
    }
    
    // Regresja liniowa dla ln(y) = c0 + c1*ln(x)
    auto linear_result = polynomial_least_squares(ln_x, ln_y, 1);
    
    if (linear_result.is_valid) {
        // Przekształcenie z powrotem: a = exp(c0), b = c1
        result.coefficients.resize(2);
        result.coefficients[0] = std::exp(linear_result.coefficients[0]); // a
        result.coefficients[1] = linear_result.coefficients[1];           // b
        result.is_valid = true;
        result = compute_approximation_statistics(x_data, y_data, result);
    }
    
    return result;
}

double Approximation::chebyshev_polynomial(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    
    // Rekurencja: T_n(x) = 2x*T_{n-1}(x) - T_{n-2}(x)
    double T_prev2 = 1.0;  // T_0(x)
    double T_prev1 = x;    // T_1(x)
    double T_current;
    
    for (int i = 2; i <= n; ++i) {
        T_current = 2.0 * x * T_prev1 - T_prev2;
        T_prev2 = T_prev1;
        T_prev1 = T_current;
    }
    
    return T_current;
}

double Approximation::map_to_chebyshev_interval(double x, double a, double b) {
    return (2.0 * x - (b + a)) / (b - a);
}

Approximation::ApproximationResult Approximation::chebyshev_approximation(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    int degree) {
    
    ApproximationResult result;
    result.method = "Chebyshev Approximation";
    result.degree = degree;
    
    if (x_data.size() != y_data.size() || x_data.empty()) {
        return result;
    }
    
    double a = *std::min_element(x_data.begin(), x_data.end());
    double b = *std::max_element(x_data.begin(), x_data.end());
    
    // Budowa funkcji bazowych Czebyszewa
    std::vector<BasisFunction> chebyshev_basis;
    for (int i = 0; i <= degree; ++i) {
        auto cheby_func = [i, a, b](double x) {
            double mapped_x = map_to_chebyshev_interval(x, a, b);
            return chebyshev_polynomial(i, mapped_x);
        };
        chebyshev_basis.emplace_back(cheby_func, "T_" + std::to_string(i));
    }
    
    // Użycie aproksymacji funkcjami bazowymi
    result = basis_function_approximation(x_data, y_data, chebyshev_basis);
    result.method = "Chebyshev Approximation";
    result.degree = degree;
    
    return result;
}

double Approximation::evaluate_approximation(const ApproximationResult& result, double x) {
    if (!result.is_valid || result.coefficients.empty()) {
        return 0.0;
    }
    
    if (result.method == "Polynomial Least Squares") {
        // Obliczenie wielomianu: P(x) = c0 + c1*x + c2*x^2 + ...
        double value = 0.0;
        double x_power = 1.0;
        for (size_t i = 0; i < result.coefficients.size(); ++i) {
            value += result.coefficients[i] * x_power;
            x_power *= x;
        }
        return value;
    }
    else if (result.method == "Exponential Approximation") {
        // f(x) = a * exp(b*x)
        return result.coefficients[0] * std::exp(result.coefficients[1] * x);
    }
    else if (result.method == "Power Approximation") {
        // f(x) = a * x^b
        if (x <= 0) return 0.0;
        return result.coefficients[0] * std::pow(x, result.coefficients[1]);
    }
    
    return 0.0; // Nieobsługiwana metoda
}

Approximation::ApproximationResult Approximation::compute_approximation_statistics(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    ApproximationResult result) {
    
    if (!result.is_valid || x_data.size() != y_data.size()) {
        return result;
    }
    
    int n = x_data.size();
    double sum_squared_errors = 0.0;
    double sum_y = 0.0;
    double sum_y_squared = 0.0;
    double max_error = 0.0;
    
    // Obliczenie średniej wartości y
    for (double y : y_data) {
        sum_y += y;
    }
    double mean_y = sum_y / n;
    
    // Obliczenie błędów i statystyk
    for (int i = 0; i < n; ++i) {
        double y_approx = evaluate_approximation(result, x_data[i]);
        double error = y_data[i] - y_approx;
        double abs_error = std::abs(error);
        
        sum_squared_errors += error * error;
        max_error = std::max(max_error, abs_error);
        
        double deviation = y_data[i] - mean_y;
        sum_y_squared += deviation * deviation;
    }
    
    // RMSE (Root Mean Square Error)
    result.rmse = std::sqrt(sum_squared_errors / n);
    
    // Maksymalny błąd bezwzględny
    result.max_error = max_error;
    
    // Współczynnik determinacji R²
    if (sum_y_squared > EPSILON) {
        result.r_squared = 1.0 - (sum_squared_errors / sum_y_squared);
    } else {
        result.r_squared = 1.0; // Wszystkie y są identyczne
    }
    
    return result;
}

int Approximation::find_optimal_polynomial_degree(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    int max_degree,
    double validation_split) {
    
    if (x_data.size() != y_data.size() || x_data.empty()) {
        return -1;
    }
    
    int n = x_data.size();
    int train_size = static_cast<int>(n * (1.0 - validation_split));
    
    if (train_size < max_degree + 1) {
        std::cerr << "Ostrzeżenie: Za mało danych treningowych dla maksymalnego stopnia" << std::endl;
        train_size = n; // Użyj wszystkich danych
    }
    
    // Podział danych na zbiór treningowy i walidacyjny
    std::vector<double> x_train(x_data.begin(), x_data.begin() + train_size);
    std::vector<double> y_train(y_data.begin(), y_data.begin() + train_size);
    std::vector<double> x_val(x_data.begin() + train_size, x_data.end());
    std::vector<double> y_val(y_data.begin() + train_size, y_data.end());
    
    double best_validation_error = std::numeric_limits<double>::max();
    int best_degree = 1;
    
    for (int degree = 1; degree <= max_degree; ++degree) {
        if (degree >= train_size) break;
        
        // Trenowanie modelu
        auto train_result = polynomial_least_squares(x_train, y_train, degree);
        
        if (!train_result.is_valid) continue;
        
        // Ewaluacja na zbiorze walidacyjnym
        double validation_error = 0.0;
        if (!x_val.empty()) {
            for (size_t i = 0; i < x_val.size(); ++i) {
                double predicted = evaluate_approximation(train_result, x_val[i]);
                double error = y_val[i] - predicted;
                validation_error += error * error;
            }
            validation_error = std::sqrt(validation_error / x_val.size());
        } else {
            // Jeśli brak danych walidacyjnych, użyj błędu treningowego
            validation_error = train_result.rmse;
        }
        
        if (validation_error < best_validation_error) {
            best_validation_error = validation_error;
            best_degree = degree;
        }
    }
    
    return best_degree;
}

std::vector<Approximation::ApproximationResult> Approximation::compare_approximation_methods(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    int max_degree) {
    
    std::vector<ApproximationResult> results;
    
    // Aproksymacja wielomianowa różnych stopni
    for (int degree = 1; degree <= max_degree; ++degree) {
        if (degree < static_cast<int>(x_data.size())) {
            auto poly_result = polynomial_least_squares(x_data, y_data, degree);
            if (poly_result.is_valid) {
                poly_result.method = "Polynomial (degree " + std::to_string(degree) + ")";
                results.push_back(poly_result);
            }
        }
    }
    
    // Aproksymacja eksponencjalna (jeśli możliwa)
    bool all_positive_y = std::all_of(y_data.begin(), y_data.end(), 
                                     [](double y) { return y > 0; });
    if (all_positive_y) {
        auto exp_result = exponential_approximation(x_data, y_data);
        if (exp_result.is_valid) {
            results.push_back(exp_result);
        }
    }
    
    // Aproksymacja potęgowa (jeśli możliwa)
    bool all_positive_x = std::all_of(x_data.begin(), x_data.end(), 
                                     [](double x) { return x > 0; });
    if (all_positive_x && all_positive_y) {
        auto power_result = power_approximation(x_data, y_data);
        if (power_result.is_valid) {
            results.push_back(power_result);
        }
    }
    
    // Aproksymacja wielomianami Czebyszewa
    int chebyshev_degree = std::min(max_degree, static_cast<int>(x_data.size()) - 1);
    auto chebyshev_result = chebyshev_approximation(x_data, y_data, chebyshev_degree);
    if (chebyshev_result.is_valid) {
        results.push_back(chebyshev_result);
    }
    
    // Sortowanie według jakości (RMSE rosnąco)
    std::sort(results.begin(), results.end(), 
              [](const ApproximationResult& a, const ApproximationResult& b) {
                  return a.rmse < b.rmse;
              });
    
    return results;
}

void Approximation::generate_approximation_plot_data(
    const ApproximationResult& result,
    double x_min, double x_max,
    int num_points,
    const std::string& filename) {
    
    if (!result.is_valid) {
        std::cerr << "Błąd: Niepoprawny wynik aproksymacji" << std::endl;
        return;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Błąd: Nie można utworzyć pliku " << filename << std::endl;
        return;
    }
    
    file << "# Dane aproksymacji - metoda: " << result.method << std::endl;
    file << "# RMSE: " << result.rmse << ", R²: " << result.r_squared << std::endl;
    file << "# x, y_approximated" << std::endl;
    
    double dx = (x_max - x_min) / (num_points - 1);
    
    for (int i = 0; i < num_points; ++i) {
        double x = x_min + i * dx;
        double y = evaluate_approximation(result, x);
        file << std::fixed << std::setprecision(8) << x << ", " << y << std::endl;
    }
    
    file.close();
    std::cout << "Dane aproksymacji zapisane do pliku: " << filename << std::endl;
}

std::vector<double> Approximation::solve_normal_equations(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& b) {
    
    // Użycie modułu układów równań liniowych
    auto solution = LinearSystems::gauss_elimination(
        const_cast<std::vector<std::vector<double>>&>(A),
        const_cast<std::vector<double>&>(b));
    
    return solution.is_valid ? solution.x : std::vector<double>();
}
Approximation::ApproximationResult Approximation::weighted_polynomial_approximation(
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<double>& weights,
    int degree) {
    
    ApproximationResult result;
    result.method = "Weighted Polynomial Approximation";
    result.degree = degree;
    
    if (x_data.size() != y_data.size() || x_data.size() != weights.size() || x_data.empty()) {
        std::cerr << "Błąd: Niezgodne rozmiary danych lub puste dane" << std::endl;
        return result;
    }
    
    if (degree < 0 || degree >= static_cast<int>(x_data.size())) {
        std::cerr << "Błąd: Niepoprawny stopień wielomianu" << std::endl;
        return result;
    }
    
    int n = x_data.size();
    int m = degree + 1; // Liczba współczynników
    
    // Budowa macierzy A układu równań normalnych z wagami: A^T * W * A * c = A^T * W * y
    std::vector<std::vector<double>> ATA(m, std::vector<double>(m, 0.0));
    std::vector<double> ATy(m, 0.0);
    
    // Wypełnienie macierzy A^T * W * A i wektora A^T * W * y
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < n; ++k) {
                ATA[i][j] += weights[k] * std::pow(x_data[k], i + j);
            }
        }
        
        for (int k = 0; k < n; ++k) {
            ATy[i] += weights[k] * y_data[k] * std::pow(x_data[k], i);
        }
    }
    
    // Rozwiązanie układu równań normalnych
    auto solution = LinearSystems::gauss_elimination(ATA, ATy);
    
    if (!solution.is_valid) {
        std::cerr << "Błąd: Nie można rozwiązać układu równań normalnych" << std::endl;
        return result;
    }
    
    result.coefficients = solution.x;
    result.is_valid = true;
    
    // Obliczenie statystyk jakości aproksymacji
    result = compute_approximation_statistics(x_data, y_data, result);
    
    return result;
}

} // namespace agh_numerical