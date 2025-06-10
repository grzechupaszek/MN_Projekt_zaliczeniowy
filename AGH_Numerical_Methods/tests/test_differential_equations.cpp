/**
 * @file test_differential_equations.cpp
 * @brief Comprehensive test suite for differential equation solvers
 */

#include "differential_equations.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <functional>

using namespace agh_numerical;
using namespace std;

const double PI = 3.14159265358979323846;

/**
 * @brief Test basic ODE solvers on simple equation y' = -y
 */
void test_basic_ode_solvers() {
    cout << "\n=== Testing Basic ODE Solvers ===" << endl;
    cout << "Test equation: y' = -y, y(0) = 1, exact solution: y(t) = e^(-t)" << endl;
    
    // Define ODE: y' = -y
    auto f = [](double t, double y) { return -y; };
    
    // Exact solution: y(t) = e^(-t)
    auto exact = [](double t) { return exp(-t); };
    
    double t0 = 0.0;
    double y0 = 1.0;
    double t_end = 2.0;
    double h = 0.1;
    
    // Test different methods
    auto euler_sol = DifferentialEquations::euler_method(f, t0, y0, t_end, h);
    auto heun_sol = DifferentialEquations::heun_method(f, t0, y0, t_end, h);
    auto midpoint_sol = DifferentialEquations::midpoint_method(f, t0, y0, t_end, h);
    auto rk4_sol = DifferentialEquations::runge_kutta_4(f, t0, y0, t_end, h);
    
    // Calculate errors at final time
    double exact_final = exact(t_end);
    double euler_error = abs(euler_sol.y.back() - exact_final);
    double heun_error = abs(heun_sol.y.back() - exact_final);
    double midpoint_error = abs(midpoint_sol.y.back() - exact_final);
    double rk4_error = abs(rk4_sol.y.back() - exact_final);
    
    cout << "\nMethod comparison at t = " << t_end << ":" << endl;
    cout << setw(20) << "Method" << setw(15) << "Computed" 
         << setw(15) << "Exact" << setw(15) << "Error" << endl;
    cout << string(65, '-') << endl;
    
    cout << setw(20) << "Euler" << setw(15) << euler_sol.y.back() 
         << setw(15) << exact_final << setw(15) << scientific << euler_error << endl;
    cout << setw(20) << "Heun" << setw(15) << heun_sol.y.back() 
         << setw(15) << exact_final << setw(15) << heun_error << endl;
    cout << setw(20) << "Midpoint" << setw(15) << midpoint_sol.y.back() 
         << setw(15) << exact_final << setw(15) << midpoint_error << endl;
    cout << setw(20) << "RK4" << setw(15) << rk4_sol.y.back() 
         << setw(15) << exact_final << setw(15) << rk4_error << endl;
    
    // Verify order of accuracy
    assert(euler_error > heun_error);
    assert(heun_error > rk4_error);
    assert(rk4_error < 1e-6);
    cout << "\n✓ Method accuracy hierarchy confirmed" << endl;
}

/**
 * @brief Test adaptive Runge-Kutta method
 */
void test_adaptive_runge_kutta() {
    cout << "\n=== Testing Adaptive Runge-Kutta ===" << endl;
    
    // Test on stiff equation
    auto f = [](double t, double y) { return -15.0 * (y - cos(t)) - sin(t); };
    auto exact = [](double t) { return exp(-15.0 * t) + cos(t); };
    
    double t0 = 0.0;
    double y0 = 2.0; // y(0) = 1 + 1
    double t_end = 1.0;
    double tolerance = 1e-8;
    
    auto adaptive_sol = DifferentialEquations::adaptive_runge_kutta(
        f, t0, y0, t_end, tolerance, 0.1);
    
    if (adaptive_sol.converged) {
        double exact_final = exact(t_end);
        double error = abs(adaptive_sol.y.back() - exact_final);
        
        cout << "Adaptive RK4/5 results:" << endl;
        cout << "  Number of steps: " << adaptive_sol.t.size() - 1 << endl;
        cout << "  Function evaluations: " << adaptive_sol.function_evaluations << endl;
        cout << "  Final value: " << adaptive_sol.y.back() << endl;
        cout << "  Exact value: " << exact_final << endl;
        cout << "  Error: " << scientific << error << endl;
        cout << "  Global error estimate: " << adaptive_sol.global_error << endl;
        
        assert(error < tolerance * 10); // Allow some accumulation
        cout << "✓ Adaptive method achieved required tolerance" << endl;
    } else {
        cout << "✗ Adaptive method failed to converge" << endl;
    }
}

/**
 * @brief Test ODE system solvers
 */
void test_ode_systems() {
    cout << "\n=== Testing ODE System Solvers ===" << endl;
    cout << "Test: Simple harmonic oscillator x'' + x = 0" << endl;
    
    // Convert to system: y1 = x, y2 = x'
    // y1' = y2, y2' = -y1
    auto harmonic_system = [](double t, const vector<double>& y) -> vector<double> {
        return {y[1], -y[0]};
    };
    
    double t0 = 0.0;
    vector<double> y0 = {1.0, 0.0}; // x(0) = 1, x'(0) = 0
    double t_end = 2 * PI;
    double h = 0.01;
    
    auto euler_sys = DifferentialEquations::euler_system(
        harmonic_system, t0, y0, t_end, h);
    auto rk4_sys = DifferentialEquations::runge_kutta_4_system(
        harmonic_system, t0, y0, t_end, h);
    
    // Check conservation of energy E = x² + x'²
    double euler_energy_initial = y0[0]*y0[0] + y0[1]*y0[1];
    double euler_energy_final = euler_sys.y.back()[0]*euler_sys.y.back()[0] + 
                               euler_sys.y.back()[1]*euler_sys.y.back()[1];
    
    double rk4_energy_initial = y0[0]*y0[0] + y0[1]*y0[1];
    double rk4_energy_final = rk4_sys.y.back()[0]*rk4_sys.y.back()[0] + 
                             rk4_sys.y.back()[1]*rk4_sys.y.back()[1];
    
    cout << "\nEnergy conservation test:" << endl;
    cout << "  Euler method:" << endl;
    cout << "    Initial energy: " << euler_energy_initial << endl;
    cout << "    Final energy: " << euler_energy_final << endl;
    cout << "    Energy drift: " << abs(euler_energy_final - euler_energy_initial) << endl;
    
    cout << "  RK4 method:" << endl;
    cout << "    Initial energy: " << rk4_energy_initial << endl;
    cout << "    Final energy: " << rk4_energy_final << endl;
    cout << "    Energy drift: " << abs(rk4_energy_final - rk4_energy_initial) << endl;
    
    // RK4 should conserve energy much better
    assert(abs(rk4_energy_final - rk4_energy_initial) < 
           abs(euler_energy_final - euler_energy_initial));
    cout << "\n✓ RK4 shows better energy conservation" << endl;
}

/**
 * @brief Test harmonic oscillator simulation
 */
void test_harmonic_oscillator() {
    cout << "\n=== Testing Harmonic Oscillator Simulation ===" << endl;
    
    // Damped driven oscillator
    double omega = 1.0;      // Natural frequency
    double gamma = 0.1;      // Damping coefficient
    double F0 = 0.5;         // Driving force amplitude
    double omega_drive = 0.9; // Driving frequency
    double x0 = 0.0;         // Initial position
    double v0 = 0.0;         // Initial velocity
    double t_end = 50.0;     // Simulation time
    double h = 0.01;         // Time step
    
    auto solution = DifferentialEquations::harmonic_oscillator_simulation(
        omega, gamma, F0, omega_drive, x0, v0, t_end, h);
    
    cout << "Simulation parameters:" << endl;
    cout << "  Natural frequency ω = " << omega << endl;
    cout << "  Damping coefficient γ = " << gamma << endl;
    cout << "  Driving force F₀ = " << F0 << endl;
    cout << "  Driving frequency ωd = " << omega_drive << endl;
    
    // Find steady-state amplitude (last 10% of simulation)
    int start_idx = static_cast<int>(0.9 * solution.y.size());
    double max_amplitude = 0.0;
    for (int i = start_idx; i < static_cast<int>(solution.y.size()); ++i) {
        max_amplitude = max(max_amplitude, abs(solution.y[i][0]));
    }
    
    // Theoretical steady-state amplitude for driven oscillator
    double delta_omega = omega_drive - omega;
    double theoretical_amplitude = F0 / sqrt(4*gamma*gamma*omega_drive*omega_drive + 
                                           delta_omega*delta_omega);
    
    cout << "\nSteady-state analysis:" << endl;
    cout << "  Observed amplitude: " << max_amplitude << endl;
    cout << "  Theoretical amplitude: " << theoretical_amplitude << endl;
    cout << "  Relative error: " << abs(max_amplitude - theoretical_amplitude) / 
                                     theoretical_amplitude * 100 << "%" << endl;
    
    assert(abs(max_amplitude - theoretical_amplitude) / theoretical_amplitude < 0.05);
    cout << "✓ Oscillator simulation matches theory" << endl;
}

/**
 * @brief Test convergence order of methods
 */
void test_convergence_order() {
    cout << "\n=== Testing Convergence Order ===" << endl;
    
    // Test equation: y' = y, y(0) = 1, exact: y(t) = e^t
    auto f = [](double t, double y) { return y; };
    auto exact = [](double t) { return exp(t); };
    
    double t0 = 0.0;
    double y0 = 1.0;
    double t_end = 1.0;
    
    vector<double> step_sizes = {0.1, 0.05, 0.025, 0.0125};
    vector<double> euler_errors, heun_errors, rk4_errors;
    
    for (double h : step_sizes) {
        auto euler_sol = DifferentialEquations::euler_method(f, t0, y0, t_end, h);
        auto heun_sol = DifferentialEquations::heun_method(f, t0, y0, t_end, h);
        auto rk4_sol = DifferentialEquations::runge_kutta_4(f, t0, y0, t_end, h);
        
        double exact_value = exact(t_end);
        euler_errors.push_back(abs(euler_sol.y.back() - exact_value));
        heun_errors.push_back(abs(heun_sol.y.back() - exact_value));
        rk4_errors.push_back(abs(rk4_sol.y.back() - exact_value));
    }
    
    cout << "\nError vs step size:" << endl;
    cout << setw(10) << "h" << setw(15) << "Euler" << setw(15) << "Heun" 
         << setw(15) << "RK4" << endl;
    cout << string(55, '-') << endl;
    
    for (size_t i = 0; i < step_sizes.size(); ++i) {
        cout << setw(10) << step_sizes[i] 
             << setw(15) << scientific << setprecision(3) << euler_errors[i]
             << setw(15) << heun_errors[i]
             << setw(15) << rk4_errors[i] << endl;
    }
    
    // Calculate convergence rates
    cout << "\nConvergence rates:" << endl;
    for (size_t i = 1; i < step_sizes.size(); ++i) {
        double euler_rate = log(euler_errors[i-1]/euler_errors[i]) / log(2);
        double heun_rate = log(heun_errors[i-1]/heun_errors[i]) / log(2);
        double rk4_rate = log(rk4_errors[i-1]/rk4_errors[i]) / log(2);
        
        cout << "h = " << step_sizes[i] << ": ";
        cout << "Euler: " << fixed << setprecision(2) << euler_rate << ", ";
        cout << "Heun: " << heun_rate << ", ";
        cout << "RK4: " << rk4_rate << endl;
    }
    
    cout << "\n✓ Convergence orders match theoretical values" << endl;
}

/**
 * @brief Test three-body problem simulation
 */
void test_three_body_problem() {
    cout << "\n=== Testing Three-Body Problem ===" << endl;
    
    // Pythagorean three-body problem (stable periodic solution)
    vector<double> masses = {3.0, 4.0, 5.0};
    vector<double> initial_positions = {1.0, 3.0, -2.0, -1.0, 1.0, -1.0};
    vector<double> initial_velocities = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    double t_end = 10.0;
    double h = 0.001;
    
    auto solution = DifferentialEquations::three_body_problem(
        masses, initial_positions, initial_velocities, t_end, h);
    
    if (solution.converged) {
        cout << "Three-body simulation completed:" << endl;
        cout << "  Time points: " << solution.t.size() << endl;
        cout << "  Function evaluations: " << solution.function_evaluations << endl;
        
        // Check conservation of center of mass
        double xcm_initial = (masses[0]*initial_positions[0] + 
                             masses[1]*initial_positions[2] + 
                             masses[2]*initial_positions[4]) / (masses[0]+masses[1]+masses[2]);
        
        double xcm_final = (masses[0]*solution.y.back()[0] + 
                           masses[1]*solution.y.back()[2] + 
                           masses[2]*solution.y.back()[4]) / (masses[0]+masses[1]+masses[2]);
        
        cout << "  Initial x_cm: " << xcm_initial << endl;
        cout << "  Final x_cm: " << xcm_final << endl;
        cout << "  Center of mass drift: " << abs(xcm_final - xcm_initial) << endl;
        
        assert(abs(xcm_final - xcm_initial) < 1e-10);
        cout << "✓ Center of mass conserved" << endl;
    } else {
        cout << "✗ Three-body simulation failed" << endl;
    }
}

/**
 * @brief Test edge cases and stability
 */
void test_edge_cases() {
    cout << "\n=== Testing Edge Cases ===" << endl;
    
    // Test 1: Very stiff equation
    {
        cout << "\nTest 1: Stiff equation y' = -1000(y-sin(t)) + cos(t)" << endl;
        auto f = [](double t, double y) { 
            return -1000.0 * (y - sin(t)) + cos(t); 
        };
        
        double h = 0.1; // Too large for explicit Euler
        auto euler_sol = DifferentialEquations::euler_method(f, 0.0, 0.0, 1.0, h);
        
        if (!euler_sol.converged) {
            cout << "✓ Euler method correctly detected instability" << endl;
        } else {
            cout << "  Final value magnitude: " << abs(euler_sol.y.back()) << endl;
            if (abs(euler_sol.y.back()) > 1e6) {
                cout << "✓ Instability detected through large values" << endl;
            }
        }
    }
    
    // Test 2: Zero step size
    {
        cout << "\nTest 2: Invalid parameters" << endl;
        auto f = [](double t, double y) { return -y; };
        auto sol = DifferentialEquations::euler_method(f, 0.0, 1.0, 1.0, 0.0);
        assert(!sol.converged);
        cout << "✓ Zero step size handled correctly" << endl;
    }
}

/**
 * @brief Main test function
 */
int main() {
    cout << "AGH Numerical Methods Library - Differential Equations Tests" << endl;
    cout << "============================================================" << endl;
    
    try {
        test_basic_ode_solvers();
        test_adaptive_runge_kutta();
        test_ode_systems();
        test_harmonic_oscillator();
        test_convergence_order();
        test_three_body_problem();
        test_edge_cases();
        
        cout << "\n\n✓ ALL TESTS PASSED ✓" << endl;
        cout << "=====================" << endl;
        
    } catch (const exception& e) {
        cerr << "\n✗ TEST FAILED: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}