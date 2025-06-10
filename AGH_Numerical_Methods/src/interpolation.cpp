/**
 * @file interpolation.cpp
 * @brief Implementacja metod interpolacji wielomianowej
 */

#include "interpolation.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <numeric>

namespace agh_numerical {

const double Interpolation::EPSILON = 1e-14;

double Interpolation::lagrange_interpolation(const std::vector<double>& nodes,
                                           const std::vector<double>& values,
                                           double x) {
    if (nodes.size() != values.size() || nodes.empty()) {
        std::cerr << "Błąd: Niepoprawne rozmiary wektorów węzłów i wartości" << std::endl;
        return 0.0;
    }
    
    double result = 0.0;
    int n = nodes.size();
    
    for (int i = 0; i < n; i++) {
        double term = values[i];
        
        // Obliczenie bazowej funkcji Lagrange'a L_i(x)
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double denominator = nodes[i] - nodes[j];
                if (std::abs(denominator) < EPSILON) {
                    std::cerr << "Ostrzeżenie: Węzły interpolacji są zbyt blisko siebie" << std::endl;
                    return 0.0;
                }
                term *= (x - nodes[j]) / denominator;
            }
        }
        
        result += term;
    }
    
    return result;
}

std::vector<std::vector<double>> Interpolation::compute_divided_differences(
    const std::vector<double>& nodes,
    const std::vector<double>& values) {
    
    int n = nodes.size();
    std::vector<std::vector<double>> table(n, std::vector<double>(n, 0.0));
    
    // Inicjalizacja pierwszej kolumny (f[x_i])
    for (int i = 0; i < n; ++i) {
        table[i][0] = values[i];
    }
    
    // Obliczenie różnic dzielonych wyższego rzędu
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            double denominator = nodes[i + j] - nodes[i];
            if (std::abs(denominator) < EPSILON) {
                std::cerr << "Błąd: Identyczne węzły interpolacji" << std::endl;
                return table;
            }
            table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / denominator;
        }
    }
    
    return table;
}

double Interpolation::newton_interpolation(const std::vector<double>& nodes,
                                         const std::vector<std::vector<double>>& divided_diffs,
                                         double x) {
    if (nodes.empty() || divided_diffs.empty()) {
        return 0.0;
    }
    
    double result = divided_diffs[0][0];
    double product = 1.0;
    
    for (size_t i = 1; i < nodes.size(); ++i) {
        product *= (x - nodes[i - 1]);
        result += divided_diffs[0][i] * product;
    }
    
    return result;
}

double Interpolation::horner_evaluation(const std::vector<double>& coefficients,
                                       double x) {
    if (coefficients.empty()) {
        return 0.0;
    }
    
    // Algorytm Hornera: P(x) = a_n + x(a_{n-1} + x(a_{n-2} + ... + x*a_1)...)
    double result = coefficients[0];
    for (size_t i = 1; i < coefficients.size(); ++i) {
        result = result * x + coefficients[i];
    }
    
    return result;
}

double Interpolation::natural_form_evaluation(const std::vector<double>& coefficients,
                                            double x) {
    if (coefficients.empty()) {
        return 0.0;
    }
    
    // Forma naturalna: P(x) = a_0 + a_1*x + a_2*x² + ... + a_n*x^n
    double result = 0.0;
    double x_power = 1.0;
    
    for (size_t i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * x_power;
        x_power *= x;
    }
    
    return result;
}

void Interpolation::select_nodes(const std::vector<double>& all_nodes,
                                const std::vector<double>& all_values,
                                int step,
                                std::vector<double>& selected_nodes,
                                std::vector<double>& selected_values) {
    selected_nodes.clear();
    selected_values.clear();
    
    if (step <= 0 || all_nodes.size() != all_values.size()) {
        std::cerr << "Błąd: Niepoprawny krok wyboru węzłów lub niezgodne rozmiary" << std::endl;
        return;
    }
    
    for (size_t i = 0; i < all_nodes.size(); i += step) {
        selected_nodes.push_back(all_nodes[i]);
        selected_values.push_back(all_values[i]);
    }
}

double Interpolation::compute_mse(const std::vector<double>& all_nodes,
                                const std::vector<double>& all_values,
                                const std::vector<double>& interp_nodes,
                                const std::vector<double>& interp_values) {
    if (all_nodes.size() != all_values.size()) {
        return -1.0; // Błąd - niezgodne rozmiary
    }
    
    double error_sum = 0.0;
    int count = 0;
    
    for (size_t i = 0; i < all_nodes.size(); i++) {
        // Sprawdź czy punkt nie jest węzłem interpolacji
        bool is_node = false;
        for (size_t j = 0; j < interp_nodes.size(); j++) {
            if (std::abs(all_nodes[i] - interp_nodes[j]) < EPSILON) {
                is_node = true;
                break;
            }
        }
        
        if (!is_node) {
            double interpolated_value = lagrange_interpolation(interp_nodes, interp_values, all_nodes[i]);
            double error = interpolated_value - all_values[i];
            error_sum += error * error;
            count++;
        }
    }
    
    return (count > 0) ? std::sqrt(error_sum / count) : 0.0;
}

int Interpolation::find_optimal_nodes_count(const std::vector<double>& nodes,
                                          const std::vector<double>& values,
                                          int min_nodes,
                                          int max_nodes) {
    if (nodes.size() != values.size() || nodes.empty()) {
        return -1;
    }
    
    int optimal_count = min_nodes;
    double best_mse = std::numeric_limits<double>::max();
    
    max_nodes = std::min(max_nodes, static_cast<int>(nodes.size()));
    
    for (int num_nodes = min_nodes; num_nodes <= max_nodes; num_nodes++) {
        // Wybierz równomiernie rozłożone węzły
        std::vector<double> selected_nodes, selected_values;
        double step = static_cast<double>(nodes.size() - 1) / (num_nodes - 1);
        
        for (int i = 0; i < num_nodes; i++) {
            int index = static_cast<int>(i * step + 0.5);
            if (index >= static_cast<int>(nodes.size())) {
                index = nodes.size() - 1;
            }
            selected_nodes.push_back(nodes[index]);
            selected_values.push_back(values[index]);
        }
        
        double mse = compute_mse(nodes, values, selected_nodes, selected_values);
        
        if (mse >= 0 && mse < best_mse) {
            best_mse = mse;
            optimal_count = num_nodes;
        }
    }
    
    return optimal_count;
}

bool Interpolation::load_data(const std::string& filename,
                            std::vector<double>& nodes,
                            std::vector<double>& values) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Nie można otworzyć pliku: " << filename << std::endl;
        return false;
    }
    
    std::string line;
    nodes.clear();
    values.clear();
    
    // Wczytaj węzły (xi:)
    if (std::getline(file, line)) {
        if (line.find("xi:") != std::string::npos) {
            std::istringstream iss(line.substr(line.find(":") + 1));
            double val;
            while (iss >> val) {
                nodes.push_back(val);
            }
        } else {
            std::cerr << "Błędny format pliku - oczekiwano 'xi:'" << std::endl;
            return false;
        }
    }
    
    // Wczytaj wartości (f(xi):)
    if (std::getline(file, line)) {
        if (line.find("f(xi):") != std::string::npos) {
            std::istringstream iss(line.substr(line.find(":") + 1));
            double val;
            while (iss >> val) {
                values.push_back(val);
            }
        } else {
            std::cerr << "Błędny format pliku - oczekiwano 'f(xi):'" << std::endl;
            return false;
        }
    }
    
    if (nodes.size() != values.size()) {
        std::cerr << "Błąd: Niezgodna liczba węzłów i wartości" << std::endl;
        return false;
    }
    
    if (nodes.empty()) {
        std::cerr << "Błąd: Brak danych w pliku" << std::endl;
        return false;
    }
    
    file.close();
    return true;
}

void Interpolation::generate_interpolation_data(const std::vector<double>& nodes,
                                              const std::vector<double>& values,
                                              const std::string& filename,
                                              int num_points,
                                              const std::string& method) {
    if (nodes.empty() || values.empty() || nodes.size() != values.size()) {
        std::cerr << "Błąd: Niepoprawne dane wejściowe" << std::endl;
        return;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Nie można utworzyć pliku: " << filename << std::endl;
        return;
    }
    
    double x_min = *std::min_element(nodes.begin(), nodes.end());
    double x_max = *std::max_element(nodes.begin(), nodes.end());
    double dx = (x_max - x_min) / (num_points - 1);
    
    file << "# Dane interpolacji - metoda: " << method << std::endl;
    file << "# x, y_interpolated" << std::endl;
    
    if (method == "newton") {
        // Użyj metody Newtona
        auto divided_diffs = compute_divided_differences(nodes, values);
        
        for (int i = 0; i < num_points; i++) {
            double x = x_min + i * dx;
            double y = newton_interpolation(nodes, divided_diffs, x);
            file << std::fixed << std::setprecision(8) << x << ", " << y << std::endl;
        }
    } else {
        // Domyślnie użyj metody Lagrange'a
        for (int i = 0; i < num_points; i++) {
            double x = x_min + i * dx;
            double y = lagrange_interpolation(nodes, values, x);
            file << std::fixed << std::setprecision(8) << x << ", " << y << std::endl;
        }
    }
    
    file.close();
    std::cout << "Dane interpolacji zapisane do pliku: " << filename << std::endl;
}

std::pair<double, double> Interpolation::compare_evaluation_methods(
    const std::vector<double>& coefficients,
    double x,
    int iterations) {
    
    if (coefficients.empty()) {
        return {0.0, 0.0};
    }
    
    // Pomiar czasu dla algorytmu Hornera
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; i++) {
        volatile double result = horner_evaluation(coefficients, x);
        (void)result; // Zapobieganie optymalizacji
    }
    auto end = std::chrono::high_resolution_clock::now();
    double horner_time = std::chrono::duration<double>(end - start).count();
    
    // Pomiar czasu dla formy naturalnej
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; i++) {
        volatile double result = natural_form_evaluation(coefficients, x);
        (void)result; // Zapobieganie optymalizacji
    }
    end = std::chrono::high_resolution_clock::now();
    double natural_time = std::chrono::duration<double>(end - start).count();
    
    return {horner_time, natural_time};
}

} // namespace agh_numerical