/**
 * @file differential_equations.cpp
 * @brief Implementacja metod rozwiązywania równań różniczkowych zwyczajnych
 */

#include "differential_equations.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

namespace agh_numerical {

const double DifferentialEquations::EPSILON = 1e-12;
const double DifferentialEquations::MIN_STEP_SIZE = 1e-10;
const double DifferentialEquations::MAX_STEP_SIZE = 1.0;

DifferentialEquations::ODESolution DifferentialEquations::euler_method(
    ODEFunction f, double t0, double y0, double t_end, double h) {
    
    ODESolution solution;
    solution.method = "Euler";
    
    if (h <= 0 || t_end <= t0) {
        std::cerr << "Błąd: Niepoprawne parametry (h > 0, t_end > t0)" << std::endl;
        solution.converged = false;
        return solution;
    }
    
    double t = t0;
    double y = y0;
    
    solution.t.push_back(t);
    solution.y.push_back(y);
    solution.function_evaluations = 0;
    
    while (t < t_end) {
        // Dostosowanie kroku na końcu przedziału
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        // Metoda Eulera: y_{n+1} = y_n + h * f(t_n, y_n)
        double f_val = f(t, y);
        solution.function_evaluations++;
        
        y = y + h * f_val;
        t = t + h;
        
        solution.t.push_back(t);
        solution.y.push_back(y);
        
        // Sprawdzenie stabilności
        if (!std::isfinite(y)) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w metodzie Eulera" << std::endl;
            solution.converged = false;
            break;
        }
    }
    
    return solution;
}

DifferentialEquations::ODESolution DifferentialEquations::heun_method(
    ODEFunction f, double t0, double y0, double t_end, double h) {
    
    ODESolution solution;
    solution.method = "Heun";
    
    if (h <= 0 || t_end <= t0) {
        std::cerr << "Błąd: Niepoprawne parametry" << std::endl;
        solution.converged = false;
        return solution;
    }
    
    double t = t0;
    double y = y0;
    
    solution.t.push_back(t);
    solution.y.push_back(y);
    solution.function_evaluations = 0;
    
    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        // Predyktor (metoda Eulera)
        double f1 = f(t, y);
        double y_tilde = y + h * f1;
        
        // Korektor (metoda Heuna)
        double f2 = f(t + h, y_tilde);
        y = y + (h / 2.0) * (f1 + f2);
        t = t + h;
        
        solution.function_evaluations += 2;
        
        solution.t.push_back(t);
        solution.y.push_back(y);
        
        if (!std::isfinite(y)) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w metodzie Heuna" << std::endl;
            solution.converged = false;
            break;
        }
    }
    
    return solution;
}

DifferentialEquations::ODESolution DifferentialEquations::midpoint_method(
    ODEFunction f, double t0, double y0, double t_end, double h) {
    
    ODESolution solution;
    solution.method = "Midpoint";
    
    if (h <= 0 || t_end <= t0) {
        std::cerr << "Błąd: Niepoprawne parametry" << std::endl;
        solution.converged = false;
        return solution;
    }
    
    double t = t0;
    double y = y0;
    
    solution.t.push_back(t);
    solution.y.push_back(y);
    solution.function_evaluations = 0;
    
    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        // Metoda punktu środkowego
        double k1 = f(t, y);
        double k2 = f(t + h/2.0, y + h * k1 / 2.0);
        
        y = y + h * k2;
        t = t + h;
        
        solution.function_evaluations += 2;
        
        solution.t.push_back(t);
        solution.y.push_back(y);
        
        if (!std::isfinite(y)) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w metodzie punktu środkowego" << std::endl;
            solution.converged = false;
            break;
        }
    }
    
    return solution;
}

DifferentialEquations::ODESolution DifferentialEquations::runge_kutta_4(
    ODEFunction f, double t0, double y0, double t_end, double h) {
    
    ODESolution solution;
    solution.method = "RK4";
    
    if (h <= 0 || t_end <= t0) {
        std::cerr << "Błąd: Niepoprawne parametry" << std::endl;
        solution.converged = false;
        return solution;
    }
    
    double t = t0;
    double y = y0;
    
    solution.t.push_back(t);
    solution.y.push_back(y);
    solution.function_evaluations = 0;
    
    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        // Klasyczna metoda RK4
        double k1 = f(t, y);
        double k2 = f(t + h/2.0, y + h * k1 / 2.0);
        double k3 = f(t + h/2.0, y + h * k2 / 2.0);
        double k4 = f(t + h, y + h * k3);
        
        y = y + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + h;
        
        solution.function_evaluations += 4;
        
        solution.t.push_back(t);
        solution.y.push_back(y);
        
        if (!std::isfinite(y)) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w metodzie RK4" << std::endl;
            solution.converged = false;
            break;
        }
    }
    
    return solution;
}

DifferentialEquations::ODESolution DifferentialEquations::adaptive_runge_kutta(
    ODEFunction f, double t0, double y0, double t_end, 
    double tolerance, double h_initial) {
    
    ODESolution solution;
    solution.method = "Adaptive RK4/5";
    
    if (tolerance <= 0 || h_initial <= 0 || t_end <= t0) {
        std::cerr << "Błąd: Niepoprawne parametry" << std::endl;
        solution.converged = false;
        return solution;
    }
    
    double t = t0;
    double y = y0;
    double h = h_initial;
    
    solution.t.push_back(t);
    solution.y.push_back(y);
    solution.function_evaluations = 0;
    solution.global_error = 0.0;
    
    const int max_iterations = 100000;
    int iterations = 0;
    
    while (t < t_end && iterations < max_iterations) {
        iterations++;
        
        // Dostosowanie kroku na końcu przedziału
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        // RK4 - metoda 4. rzędu
        double k1 = f(t, y);
        double k2 = f(t + h/2.0, y + h * k1 / 2.0);
        double k3 = f(t + h/2.0, y + h * k2 / 2.0);
        double k4 = f(t + h, y + h * k3);
        
        double y4 = y + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
        
        // RK5 - metoda 5. rzędu (uproszczona)
        double k5 = f(t + h, y4);
        double y5 = y + (h / 90.0) * (7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5);
        
        solution.function_evaluations += 5;
        
        // Oszacowanie błędu lokalnego
        double error = std::abs(y5 - y4);
        
        if (error <= tolerance || h <= MIN_STEP_SIZE) {
            // Akceptuj krok
            y = y4; // Używamy rozwiązania RK4
            t = t + h;
            
            solution.t.push_back(t);
            solution.y.push_back(y);
            solution.error.push_back(error);
            solution.global_error += error;
            
            if (!std::isfinite(y)) {
                std::cerr << "Ostrzeżenie: Niestabilność numeryczna" << std::endl;
                solution.converged = false;
                break;
            }
        }
        
        // Dostosowanie kroku na następną iterację
        if (error > EPSILON) {
            double optimal_h = compute_optimal_step_size(h, error, tolerance, 4);
            h = std::max(MIN_STEP_SIZE, std::min(MAX_STEP_SIZE, optimal_h));
        }
    }
    
    if (iterations >= max_iterations) {
        std::cerr << "Ostrzeżenie: Osiągnięto maksymalną liczbę iteracji" << std::endl;
        solution.converged = false;
    }
    
    return solution;
}

DifferentialEquations::ODESystemSolution DifferentialEquations::euler_system(
    ODESystemFunction f, double t0, const std::vector<double>& y0,
    double t_end, double h) {
    
    ODESystemSolution solution;
    solution.method = "Euler System";
    
    if (h <= 0 || t_end <= t0 || y0.empty()) {
        std::cerr << "Błąd: Niepoprawne parametry układu" << std::endl;
        solution.converged = false;
        return solution;
    }
    
    double t = t0;
    std::vector<double> y = y0;
    int n = y0.size();
    
    solution.t.push_back(t);
    solution.y.push_back(y);
    solution.function_evaluations = 0;
    
    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        // Metoda Eulera dla układu
        std::vector<double> dy_dt = f(t, y);
        solution.function_evaluations++;
        
        for (int i = 0; i < n; ++i) {
            y[i] = y[i] + h * dy_dt[i];
        }
        t = t + h;
        
        solution.t.push_back(t);
        solution.y.push_back(y);
        
        // Sprawdzenie stabilności
        bool stable = true;
        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(y[i])) {
                stable = false;
                break;
            }
        }
        
        if (!stable) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w układzie" << std::endl;
            solution.converged = false;
            break;
        }
    }
    
    return solution;
}

DifferentialEquations::ODESystemSolution DifferentialEquations::runge_kutta_4_system(
    ODESystemFunction f, double t0, const std::vector<double>& y0,
    double t_end, double h) {
    
    ODESystemSolution solution;
    solution.method = "RK4 System";
    
    if (h <= 0 || t_end <= t0 || y0.empty()) {
        std::cerr << "Błąd: Niepoprawne parametry układu" << std::endl;
        solution.converged = false;
        return solution;
    }
    
    double t = t0;
    std::vector<double> y = y0;
    int n = y0.size();
    
    solution.t.push_back(t);
    solution.y.push_back(y);
    solution.function_evaluations = 0;
    
    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }
        
        // RK4 dla układu równań
        std::vector<double> k1 = f(t, y);
        
        std::vector<double> y_temp(n);
        for (int i = 0; i < n; ++i) {
            y_temp[i] = y[i] + h * k1[i] / 2.0;
        }
        std::vector<double> k2 = f(t + h/2.0, y_temp);
        
        for (int i = 0; i < n; ++i) {
            y_temp[i] = y[i] + h * k2[i] / 2.0;
        }
        std::vector<double> k3 = f(t + h/2.0, y_temp);
        
        for (int i = 0; i < n; ++i) {
            y_temp[i] = y[i] + h * k3[i];
        }
        std::vector<double> k4 = f(t + h, y_temp);
        
        solution.function_evaluations += 4;
        
        // Aktualizacja rozwiązania
        for (int i = 0; i < n; ++i) {
            y[i] = y[i] + (h / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        }
        t = t + h;
        
        solution.t.push_back(t);
        solution.y.push_back(y);
        
        // Sprawdzenie stabilności
        bool stable = true;
        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(y[i])) {
                stable = false;
                break;
            }
        }
        
        if (!stable) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w układzie RK4" << std::endl;
            solution.converged = false;
            break;
        }
    }
    
    return solution;
}

DifferentialEquations::ODESystemSolution DifferentialEquations::harmonic_oscillator_simulation(
    double omega, double gamma, double F0, double omega_drive,
    double x0, double v0, double t_end, double h) {
    
    // Układ równań: ẍ + 2γẋ + ω²x = F₀cos(ωd·t)
    // Przekształcenie do układu 1. rzędu: y₁ = x, y₂ = ẋ
    // y₁' = y₂
    // y₂' = -2γy₂ - ω²y₁ + F₀cos(ωd·t)
    
    auto oscillator_system = [omega, gamma, F0, omega_drive](
        double t, const std::vector<double>& y) -> std::vector<double> {
        
        std::vector<double> dy_dt(2);
        dy_dt[0] = y[1]; // dx/dt = v
        dy_dt[1] = -2.0 * gamma * y[1] - omega * omega * y[0] + 
                   F0 * std::cos(omega_drive * t); // dv/dt
        
        return dy_dt;
    };
    
    std::vector<double> initial_conditions = {x0, v0};
    
    return runge_kutta_4_system(oscillator_system, 0.0, initial_conditions, t_end, h);
}

DifferentialEquations::ODESystemSolution DifferentialEquations::three_body_problem(
    const std::vector<double>& masses,
    const std::vector<double>& initial_positions,
    const std::vector<double>& initial_velocities,
    double t_end, double h) {
    
    if (masses.size() != 3 || initial_positions.size() != 6 || initial_velocities.size() != 6) {
        ODESystemSolution solution;
        solution.converged = false;
        std::cerr << "Błąd: Niepoprawne wymiary danych dla problemu trzech ciał" << std::endl;
        return solution;
    }
    
    // Stan układu: [x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3]
    auto three_body_system = [masses](double t, const std::vector<double>& state) -> std::vector<double> {
        (void)t; // Nieużywane - układ autonomiczny
        
        std::vector<double> derivative(12);
        
        // Pozycje i prędkości
        double x1 = state[0], y1 = state[1];
        double x2 = state[6], y2 = state[7];
        double x3 = state[12], y3 = state[13];
        
        double vx1 = state[6], vy1 = state[7];
        double vx2 = state[8], vy2 = state[9];
        double vx3 = state[10], vy3 = state[11];
        
        // Poprawka indeksów
        x2 = state[2]; y2 = state[3];
        x3 = state[4]; y3 = state[5];
        
        // Pochodne pozycji = prędkości
        derivative[0] = vx1; derivative[1] = vy1;
        derivative[2] = vx2; derivative[3] = vy2;
        derivative[4] = vx3; derivative[5] = vy3;
        
        // Obliczenie odległości
        double r12 = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        double r13 = std::sqrt((x3-x1)*(x3-x1) + (y3-y1)*(y3-y1));
        double r23 = std::sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2));
        
        // Unikanie singularności
        const double eps = 1e-6;
        r12 = std::max(r12, eps);
        r13 = std::max(r13, eps);
        r23 = std::max(r23, eps);
        
        double G = 1.0; // Stała grawitacji
        
        // Przyspieszenia ciała 1
        derivative[6] = G * masses[1] * (x2-x1) / (r12*r12*r12) + 
                       G * masses[2] * (x3-x1) / (r13*r13*r13);
        derivative[7] = G * masses[1] * (y2-y1) / (r12*r12*r12) + 
                       G * masses[2] * (y3-y1) / (r13*r13*r13);
        
        // Przyspieszenia ciała 2
        derivative[8] = G * masses[0] * (x1-x2) / (r12*r12*r12) + 
                       G * masses[2] * (x3-x2) / (r23*r23*r23);
        derivative[9] = G * masses[0] * (y1-y2) / (r12*r12*r12) + 
                       G * masses[2] * (y3-y2) / (r23*r23*r23);
        
        // Przyspieszenia ciała 3
        derivative[10] = G * masses[0] * (x1-x3) / (r13*r13*r13) + 
                        G * masses[1] * (x2-x3) / (r23*r23*r23);
        derivative[11] = G * masses[0] * (y1-y3) / (r13*r13*r13) + 
                        G * masses[1] * (y2-y3) / (r23*r23*r23);
        
        return derivative;
    };
    
    // Stan początkowy
    std::vector<double> initial_state(12);
    for (int i = 0; i < 6; ++i) {
        initial_state[i] = initial_positions[i];
        initial_state[i + 6] = initial_velocities[i];
    }
    
    return runge_kutta_4_system(three_body_system, 0.0, initial_state, t_end, h);
}

void DifferentialEquations::save_solution_to_csv(
    const ODESolution& solution, const std::string& filename, bool include_error) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Błąd: Nie można utworzyć pliku " << filename << std::endl;
        return;
    }
    
    file << "# Rozwiązanie równania różniczkowego - metoda: " << solution.method << std::endl;
    file << "# Wywołania funkcji: " << solution.function_evaluations << std::endl;
    file << "# Zbieżność: " << (solution.converged ? "TAK" : "NIE") << std::endl;
    
    if (include_error && !solution.error.empty()) {
        file << "t,y,error" << std::endl;
        for (size_t i = 0; i < solution.t.size(); ++i) {
            file << std::fixed << std::setprecision(8) << solution.t[i] << ","
                 << solution.y[i] << ",";
            if (i < solution.error.size()) {
                file << std::scientific << solution.error[i];
            } else {
                file << "0";
            }
            file << std::endl;
        }
    } else {
        file << "t,y" << std::endl;
        for (size_t i = 0; i < solution.t.size(); ++i) {
            file << std::fixed << std::setprecision(8) << solution.t[i] << ","
                 << solution.y[i] << std::endl;
        }
    }
    
    file.close();
    std::cout << "Rozwiązanie zapisane do pliku: " << filename << std::endl;
}

void DifferentialEquations::save_system_solution_to_csv(
    const ODESystemSolution& solution, const std::string& filename,
    const std::vector<std::string>& variable_names) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Błąd: Nie można utworzyć pliku " << filename << std::endl;
        return;
    }
    
    file << "# Rozwiązanie układu równań różniczkowych - metoda: " << solution.method << std::endl;
    file << "# Wywołania funkcji: " << solution.function_evaluations << std::endl;
    file << "# Zbieżność: " << (solution.converged ? "TAK" : "NIE") << std::endl;
    
    // Nagłówek
    file << "t";
    if (!variable_names.empty()) {
        for (const auto& name : variable_names) {
            file << "," << name;
        }
    } else {
        if (!solution.y.empty()) {
            for (size_t j = 0; j < solution.y[0].size(); ++j) {
                file << ",y" << j;
            }
        }
    }
    file << std::endl;
    
    // Dane
    for (size_t i = 0; i < solution.t.size(); ++i) {
        file << std::fixed << std::setprecision(8) << solution.t[i];
        if (i < solution.y.size()) {
            for (size_t j = 0; j < solution.y[i].size(); ++j) {
                file << "," << solution.y[i][j];
            }
        }
        file << std::endl;
    }
    
    file.close();
    std::cout << "Rozwiązanie układu zapisane do pliku: " << filename << std::endl;
}

void DifferentialEquations::convergence_analysis(
    ODEFunction f, std::function<double(double)> analytical_solution,
    double t0, double y0, double t_end, const std::string& filename) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Błąd: Nie można utworzyć pliku " << filename << std::endl;
        return;
    }
    
    file << "# Analiza zbieżności metod rozwiązywania równań różniczkowych" << std::endl;
    file << "h,Euler_error,Heun_error,RK4_error,Euler_rate,Heun_rate,RK4_rate" << std::endl;
    
    std::vector<double> step_sizes = {0.1, 0.05, 0.025, 0.0125, 0.00625};
    std::vector<double> euler_errors, heun_errors, rk4_errors;
    
    double exact_final = analytical_solution(t_end);
    
    for (double h : step_sizes) {
        auto euler_sol = euler_method(f, t0, y0, t_end, h);
        auto heun_sol = heun_method(f, t0, y0, t_end, h);
        auto rk4_sol = runge_kutta_4(f, t0, y0, t_end, h);
        
        double euler_error = std::abs(euler_sol.y.back() - exact_final);
        double heun_error = std::abs(heun_sol.y.back() - exact_final);
        double rk4_error = std::abs(rk4_sol.y.back() - exact_final);
        
        euler_errors.push_back(euler_error);
        heun_errors.push_back(heun_error);
        rk4_errors.push_back(rk4_error);
        
        file << h << "," << std::scientific << euler_error << "," 
             << heun_error << "," << rk4_error;
        
        // Obliczenie tempa zbieżności
        if (euler_errors.size() > 1) {
            double euler_rate = std::log(euler_errors[euler_errors.size()-2] / euler_error) / std::log(2.0);
            double heun_rate = std::log(heun_errors[heun_errors.size()-2] / heun_error) / std::log(2.0);
            double rk4_rate = std::log(rk4_errors[rk4_errors.size()-2] / rk4_error) / std::log(2.0);
            
            file << "," << std::fixed << std::setprecision(2) 
                 << euler_rate << "," << heun_rate << "," << rk4_rate;
        } else {
            file << ",-,-,-";
        }
        file << std::endl;
    }
    
    file.close();
    std::cout << "Analiza zbieżności zapisana do pliku: " << filename << std::endl;
}

double DifferentialEquations::compute_optimal_step_size(
    double h_current, double error, double tolerance, int order) {
    
    if (error <= EPSILON) {
        return h_current * 2.0; // Zwiększ krok jeśli błąd bardzo mały
    }
    
    // Wzór na optymalny krok: h_opt = h_current * (tolerance/error)^(1/(order+1))
    double factor = std::pow(tolerance / error, 1.0 / (order + 1));
    
    // Ogranicz zmiany kroku (czynnik bezpieczeństwa 0.9)
    factor = 0.9 * std::min(2.0, std::max(0.5, factor));
    
    return h_current * factor;
}

} // namespace agh_numerical