/**
 * @file example_differential_equations.cpp
 * @brief Przykład demonstracyjny użycia metod rozwiązywania równań różniczkowych
 * @author Student Inżynierii Obliczeniowej AGH
 * 
 * @details Ten przykład pokazuje praktyczne zastosowanie różnych metod
 *          rozwiązywania równań różniczkowych zwyczajnych dostępnych 
 *          w bibliotece AGH Numerical Methods. Demonstruje metody Eulera,
 *          Heuna, Runge-Kutta oraz ich zastosowania w problemach fizycznych.
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
 * @brief Demonstracja podstawowych metod rozwiązywania równań różniczkowych
 */
void demo_basic_ode_methods() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "DEMONSTRACJA PODSTAWOWYCH METOD ROZWIĄZYWANIA RÓR" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Równanie testowe: y' = -2y + t, y(0) = 1
    // Rozwiązanie analityczne: y(t) = (3*exp(-2t) + 2t - 1)/2
    auto differential_eq = [](double t, double y) {
        return -2.0 * y + t;
    };
    
    auto analytical_solution = [](double t) {
        return (3.0 * std::exp(-2.0 * t) + 2.0 * t - 1.0) / 2.0;
    };
    
    double t0 = 0.0, y0 = 1.0, t_end = 2.0;
    double h = 0.1;
    
    std::cout << "Równanie różniczkowe: y' = -2y + t" << std::endl;
    std::cout << "Warunek początkowy: y(0) = 1" << std::endl;
    std::cout << "Rozwiązanie analityczne: y(t) = (3*exp(-2t) + 2t - 1)/2" << std::endl;
    std::cout << "Przedział całkowania: [" << t0 << ", " << t_end << "]" << std::endl;
    std::cout << "Krok całkowania: h = " << h << std::endl;
    
    // Porównanie różnych metod
    auto euler_solution = DifferentialEquations::euler_method(differential_eq, t0, y0, t_end, h);
    auto heun_solution = DifferentialEquations::heun_method(differential_eq, t0, y0, t_end, h);
    auto midpoint_solution = DifferentialEquations::midpoint_method(differential_eq, t0, y0, t_end, h);
    auto rk4_solution = DifferentialEquations::runge_kutta_4(differential_eq, t0, y0, t_end, h);
    
    std::cout << "\n" << std::string(80, '-') << std::endl;
    std::cout << std::setw(8) << "t" << std::setw(12) << "Analityczne" 
              << std::setw(12) << "Euler" << std::setw(12) << "Heun" 
              << std::setw(12) << "Punkt śr." << std::setw(12) << "RK4" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    // Porównanie wyników w wybranych punktach
    for (size_t i = 0; i < euler_solution.t.size(); i += 5) { // Co 5-ty punkt
        double t = euler_solution.t[i];
        double exact = analytical_solution(t);
        
        std::cout << std::fixed << std::setprecision(2) << std::setw(8) << t
                  << std::setprecision(6) << std::setw(12) << exact
                  << std::setw(12) << euler_solution.y[i]
                  << std::setw(12) << heun_solution.y[i]
                  << std::setw(12) << midpoint_solution.y[i]
                  << std::setw(12) << rk4_solution.y[i] << std::endl;
    }
    
    // Analiza błędów końcowych
    double exact_final = analytical_solution(t_end);
    
    std::cout << "\nBłędy końcowe w t = " << t_end << ":" << std::endl;
    std::cout << "  Euler:        " << std::scientific << std::setprecision(3) 
              << std::abs(euler_solution.y.back() - exact_final) << std::endl;
    std::cout << "  Heun:         " << std::abs(heun_solution.y.back() - exact_final) << std::endl;
    std::cout << "  Punkt środka: " << std::abs(midpoint_solution.y.back() - exact_final) << std::endl;
    std::cout << "  RK4:          " << std::abs(rk4_solution.y.back() - exact_final) << std::endl;
    
    // Zapisanie wyników do pliku
    DifferentialEquations::save_solution_to_csv(rk4_solution, "ode_basic_example.csv");
    std::cout << "\nWyniki RK4 zapisane do pliku: ode_basic_example.csv" << std::endl;
}

/**
 * @brief Demonstracja adaptacyjnej metody Runge-Kutta
 */
void demo_adaptive_runge_kutta() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "DEMONSTRACJA ADAPTACYJNEJ METODY RUNGE-KUTTA" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Równanie o zmiennej dynamice: y' = -10y + 10sin(10t), y(0) = 0
    auto stiff_equation = [](double t, double y) {
        return -10.0 * y + 10.0 * std::sin(10.0 * t);
    };
    
    double t0 = 0.0, y0 = 0.0, t_end = 2.0;
    double tolerance = 1e-6;
    
    std::cout << "Równanie różniczkowe: y' = -10y + 10sin(10t)" << std::endl;
    std::cout << "Warunek początkowy: y(0) = 0" << std::endl;
    std::cout << "Charakterystyka: szybko zmieniająca się dynamika" << std::endl;
    std::cout << "Żądana tolerancja: " << std::scientific << tolerance << std::endl;
    
    // Metoda ze stałym krokiem
    auto fixed_step_solution = DifferentialEquations::runge_kutta_4(stiff_equation, t0, y0, t_end, 0.01);
    
    // Metoda adaptacyjna
    auto adaptive_solution = DifferentialEquations::adaptive_runge_kutta(
        stiff_equation, t0, y0, t_end, tolerance, 0.1);
    
    std::cout << "\nPorównanie metod:" << std::endl;
    std::cout << "1. RK4 ze stałym krokiem h = 0.01:" << std::endl;
    std::cout << "   Liczba kroków: " << fixed_step_solution.t.size() << std::endl;
    std::cout << "   Wywołania funkcji: " << fixed_step_solution.function_evaluations << std::endl;
    std::cout << "   Wartość końcowa: " << std::fixed << std::setprecision(8) 
              << fixed_step_solution.y.back() << std::endl;
    
    std::cout << "\n2. Adaptacyjna RK4/5:" << std::endl;
    std::cout << "   Liczba kroków: " << adaptive_solution.t.size() << std::endl;
    std::cout << "   Wywołania funkcji: " << adaptive_solution.function_evaluations << std::endl;
    std::cout << "   Wartość końcowa: " << adaptive_solution.y.back() << std::endl;
    std::cout << "   Oszacowany błąd globalny: " << std::scientific 
              << adaptive_solution.global_error << std::endl;
    std::cout << "   Zbieżność: " << (adaptive_solution.converged ? "TAK" : "NIE") << std::endl;
    
    // Efektywność
    double efficiency_gain = static_cast<double>(fixed_step_solution.function_evaluations) / 
                           adaptive_solution.function_evaluations;
    
    std::cout << "\n3. Analiza efektywności:" << std::endl;
    std::cout << "   Redukcja wywołań funkcji: " << std::fixed << std::setprecision(1)
              << (1.0 - 1.0/efficiency_gain) * 100 << "%" << std::endl;
    std::cout << "   Współczynnik efektywności: " << std::setprecision(2) 
              << efficiency_gain << "x" << std::endl;
    
    // Zapisanie wyników
    DifferentialEquations::save_solution_to_csv(adaptive_solution, "ode_adaptive_example.csv", true);
    std::cout << "\nWyniki adaptacyjne zapisane do pliku: ode_adaptive_example.csv" << std::endl;
}

/**
 * @brief Symulacja oscylatora harmonicznego
 */
void demo_harmonic_oscillator() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "SYMULACJA OSCYLATORA HARMONICZNEGO" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Parametry oscylatora: ẍ + 2γẋ + ω²x = F₀cos(ωdt)
    double omega = 2.0;        // Częstość własna [rad/s]
    double gamma = 0.1;        // Współczynnik tłumienia [1/s]
    double F0 = 1.0;           // Amplituda siły wymuszającej [N]
    double omega_drive = 1.8;  // Częstość wymuszająca [rad/s]
    double x0 = 1.0;           // Początkowe położenie [m]
    double v0 = 0.0;           // Początkowa prędkość [m/s]
    double t_end = 20.0;       // Czas symulacji [s]
    double h = 0.01;           // Krok czasowy [s]
    
    std::cout << "Równanie oscylatora: ẍ + 2γẋ + ω²x = F₀cos(ωd·t)" << std::endl;
    std::cout << "Parametry:" << std::endl;
    std::cout << "  ω = " << omega << " rad/s (częstość własna)" << std::endl;
    std::cout << "  γ = " << gamma << " 1/s (tłumienie)" << std::endl;
    std::cout << "  F₀ = " << F0 << " N (amplituda siły)" << std::endl;
    std::cout << "  ωd = " << omega_drive << " rad/s (częstość wymuszająca)" << std::endl;
    std::cout << "Warunki początkowe: x(0) = " << x0 << " m, v(0) = " << v0 << " m/s" << std::endl;
    
    // Symulacja oscylatora
    auto oscillator_solution = DifferentialEquations::harmonic_oscillator_simulation(
        omega, gamma, F0, omega_drive, x0, v0, t_end, h);
    
    std::cout << "\nWyniki symulacji:" << std::endl;
    std::cout << "  Liczba kroków czasowych: " << oscillator_solution.t.size() << std::endl;
    std::cout << "  Wywołania funkcji: " << oscillator_solution.function_evaluations << std::endl;
    
    // Analiza amplitudy w stanie ustalonym (ostatnie 25% czasu)
    size_t steady_state_start = 3 * oscillator_solution.t.size() / 4;
    double max_amplitude = 0.0;
    double min_amplitude = 0.0;
    
    for (size_t i = steady_state_start; i < oscillator_solution.t.size(); ++i) {
        double x = oscillator_solution.y[i][0]; // Pozycja
        max_amplitude = std::max(max_amplitude, x);
        min_amplitude = std::min(min_amplitude, x);
    }
    
    double amplitude = (max_amplitude - min_amplitude) / 2.0;
    
    std::cout << "  Amplituda w stanie ustalonym: " << std::fixed << std::setprecision(4) 
              << amplitude << " m" << std::endl;
    
    // Teoretyczna amplituda rezonansu
    double omega_ratio = omega_drive / omega;
    double theoretical_amplitude = F0 / (omega * omega * std::sqrt(
        std::pow(1 - omega_ratio * omega_ratio, 2) + 
        std::pow(2 * gamma * omega_ratio / omega, 2)));
    
    std::cout << "  Teoretyczna amplituda: " << theoretical_amplitude << " m" << std::endl;
    std::cout << "  Błąd względny: " << std::setprecision(2) 
              << std::abs(amplitude - theoretical_amplitude) / theoretical_amplitude * 100 
              << "%" << std::endl;
    
    // Zapisanie wyników
    DifferentialEquations::save_system_solution_to_csv(
        oscillator_solution, "oscillator_simulation.csv", {"pozycja", "predkosc"});
    std::cout << "\nWyniki symulacji zapisane do pliku: oscillator_simulation.csv" << std::endl;
    
    // Analiza energii (dla przypadku bez tłumienia i wymuszania)
    if (gamma == 0.0 && F0 == 0.0) {
        std::cout << "\nAnaliza zachowania energii (oscylator bez tłumienia):" << std::endl;
        
        double initial_energy = 0.5 * omega * omega * x0 * x0 + 0.5 * v0 * v0;
        double final_x = oscillator_solution.y.back()[0];
        double final_v = oscillator_solution.y.back()[1];
        double final_energy = 0.5 * omega * omega * final_x * final_x + 0.5 * final_v * final_v;
        
        std::cout << "  Energia początkowa: " << std::setprecision(6) << initial_energy << " J" << std::endl;
        std::cout << "  Energia końcowa: " << final_energy << " J" << std::endl;
        std::cout << "  Błąd zachowania energii: " << std::scientific 
                  << std::abs(final_energy - initial_energy) << " J" << std::endl;
    }
}

/**
 * @brief Symulacja problemu trzech ciał (ograniczona wersja)
 */
void demo_three_body_problem() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "SYMULACJA PROBLEMU TRZECH CIAŁ" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Konfiguracja trzech ciał w układzie trójkątnym
    std::vector<double> masses = {1.0, 1.0, 1.0};  // Masy [jednostki arbitralne]
    
    // Pozycje początkowe (wierzchołki trójkąta równobocznego)
    double L = 1.0; // Długość boku trójkąta
    std::vector<double> initial_positions = {
        0.0, 0.0,                           // Ciało 1
        L, 0.0,                             // Ciało 2  
        L/2, L * std::sqrt(3.0)/2           // Ciało 3
    };
    
    // Prędkości początkowe (ruch orbitalny)
    double v_orbital = 0.5;
    std::vector<double> initial_velocities = {
        0.0, v_orbital,                     // Ciało 1
        -v_orbital * std::sqrt(3.0)/2, -v_orbital/2,  // Ciało 2
        v_orbital * std::sqrt(3.0)/2, -v_orbital/2    // Ciało 3
    };
    
    double t_end = 10.0;   // Czas symulacji
    double h = 0.001;      // Mały krok dla stabilności
    
    std::cout << "Konfiguracja: trzy ciała o równych masach" << std::endl;
    std::cout << "Układ początkowy: trójkąt równoboczny z ruchem orbitalnym" << std::endl;
    std::cout << "Czas symulacji: " << t_end << " jednostek czasu" << std::endl;
    std::cout << "Krok czasowy: " << h << std::endl;
    
    std::cout << "\nRozpoczynanie symulacji..." << std::endl;
    
    auto three_body_solution = DifferentialEquations::three_body_problem(
        masses, initial_positions, initial_velocities, t_end, h);
    
    std::cout << "Symulacja zakończona." << std::endl;
    std::cout << "Liczba kroków czasowych: " << three_body_solution.t.size() << std::endl;
    std::cout << "Wywołania funkcji: " << three_body_solution.function_evaluations << std::endl;
    std::cout << "Zbieżność: " << (three_body_solution.converged ? "TAK" : "NIE") << std::endl;
    
    // Analiza zachowania środka masy
    std::cout << "\nAnaliza zachowania środka masy:" << std::endl;
    
    // Środek masy na początku
    double total_mass = masses[0] + masses[1] + masses[2];
    double cm_x_initial = (masses[0] * initial_positions[0] + 
                          masses[1] * initial_positions[2] + 
                          masses[2] * initial_positions[4]) / total_mass;
    double cm_y_initial = (masses[0] * initial_positions[1] + 
                          masses[1] * initial_positions[3] + 
                          masses[2] * initial_positions[5]) / total_mass;
    
    // Środek masy na końcu
    auto& final_state = three_body_solution.y.back();
    double cm_x_final = (masses[0] * final_state[0] + 
                        masses[1] * final_state[6] + 
                        masses[2] * final_state[12]) / total_mass;
    double cm_y_final = (masses[0] * final_state[1] + 
                        masses[1] * final_state[7] + 
                        masses[2] * final_state[13]) / total_mass;
    
    std::cout << "  Początkowy środek masy: (" << std::fixed << std::setprecision(6) 
              << cm_x_initial << ", " << cm_y_initial << ")" << std::endl;
    std::cout << "  Końcowy środek masy: (" << cm_x_final << ", " << cm_y_final << ")" << std::endl;
    std::cout << "  Przesunięcie środka masy: " << std::scientific 
              << std::sqrt(std::pow(cm_x_final - cm_x_initial, 2) + 
                          std::pow(cm_y_final - cm_y_initial, 2)) << std::endl;
    
    // Zapisanie wyników
    std::vector<std::string> variable_names = {
        "x1", "y1", "vx1", "vy1",
        "x2", "y2", "vx2", "vy2", 
        "x3", "y3", "vx3", "vy3"
    };
    
    DifferentialEquations::save_system_solution_to_csv(
        three_body_solution, "three_body_simulation.csv", variable_names);
    std::cout << "\nWyniki symulacji zapisane do pliku: three_body_simulation.csv" << std::endl;
}

/**
 * @brief Analiza zbieżności metod numerycznych
 */
void demo_convergence_analysis() {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "ANALIZA ZBIEŻNOŚCI METOD NUMERYCZNYCH" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Równanie testowe: y' = -y, y(0) = 1
    // Rozwiązanie analityczne: y(t) = exp(-t)
    auto simple_equation = [](double t, double y) {
        (void)t; // Nie używane w tym równaniu
        return -y;
    };
    
    auto analytical = [](double t) {
        return std::exp(-t);
    };
    
    double t0 = 0.0, y0 = 1.0, t_end = 1.0;
    
    std::cout << "Równanie testowe: y' = -y, y(0) = 1" << std::endl;
    std::cout << "Rozwiązanie analityczne: y(t) = exp(-t)" << std::endl;
    std::cout << "Test zbieżności w punkcie t = " << t_end << std::endl;
    std::cout << "Dokładna wartość: y(" << t_end << ") = " << std::fixed << std::setprecision(8) 
              << analytical(t_end) << std::endl;
    
    // Test dla różnych kroków czasowych
    std::vector<double> step_sizes = {0.1, 0.05, 0.025, 0.0125, 0.00625};
    
    std::cout << "\n" << std::string(70, '-') << std::endl;
    std::cout << std::setw(12) << "Krok h" << std::setw(15) << "Euler" 
              << std::setw(15) << "Heun" << std::setw(15) << "RK4" 
              << std::setw(15) << "Tempo RK4" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    double exact_value = analytical(t_end);
    double previous_rk4_error = 0.0;
    
    for (size_t i = 0; i < step_sizes.size(); ++i) {
        double h = step_sizes[i];
        
        auto euler_sol = DifferentialEquations::euler_method(simple_equation, t0, y0, t_end, h);
        auto heun_sol = DifferentialEquations::heun_method(simple_equation, t0, y0, t_end, h);
        auto rk4_sol = DifferentialEquations::runge_kutta_4(simple_equation, t0, y0, t_end, h);
        
        double euler_error = std::abs(euler_sol.y.back() - exact_value);
        double heun_error = std::abs(heun_sol.y.back() - exact_value);
        double rk4_error = std::abs(rk4_sol.y.back() - exact_value);
        
        double convergence_rate = 0.0;
        if (i > 0 && previous_rk4_error > 0) {
            convergence_rate = std::log(previous_rk4_error / rk4_error) / std::log(2.0);
        }
        
        std::cout << std::setw(12) << std::fixed << std::setprecision(5) << h
                  << std::setw(15) << std::scientific << std::setprecision(3) << euler_error
                  << std::setw(15) << heun_error
                  << std::setw(15) << rk4_error
                  << std::setw(15) << std::fixed << std::setprecision(2);
        
        if (i > 0) {
            std::cout << convergence_rate;
        } else {
            std::cout << "-";
        }
        std::cout << std::endl;
        
        previous_rk4_error = rk4_error;
    }
    
    std::cout << "\nWnioski z analizy zbieżności:" << std::endl;
    std::cout << "- Metoda Eulera: zbieżność O(h) - błąd maleje liniowo" << std::endl;
    std::cout << "- Metoda Heuna: zbieżność O(h²) - błąd maleje kwadratowo" << std::endl;
    std::cout << "- Metoda RK4: zbieżność O(h⁴) - tempo ~4 dla małych h" << std::endl;
    
    // Generacja szczegółowych danych zbieżności
    const std::string filename = "ode_convergence_analysis.csv";
    DifferentialEquations::convergence_analysis(simple_equation, analytical, t0, y0, t_end, filename);
    std::cout << "\nSzczegółowa analiza zbieżności zapisana do pliku: " << filename << std::endl;
}

/**
 * @brief Funkcja główna - demonstracja wszystkich funkcjonalności
 */
int main() {
    std::cout << "AGH NUMERICAL METHODS LIBRARY - PRZYKŁAD RÓWNAŃ RÓŻNICZKOWYCH" << std::endl;
    std::cout << "Biblioteka metod numerycznych dla inżynierii obliczeniowej" << std::endl;
    std::cout << "Autor: Student Inżynierii Obliczeniowej AGH" << std::endl;
    
    try {
        // Demonstracje poszczególnych funkcjonalności
        demo_basic_ode_methods();
        demo_adaptive_runge_kutta();
        demo_harmonic_oscillator();
        demo_three_body_problem();
        demo_convergence_analysis();
        
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "DEMONSTRACJA ZAKOŃCZONA POMYŚLNIE" << std::endl;
        std::cout << std::string(70, '=') << std::endl;
        
        std::cout << "\nPodsumowanie funkcjonalności:" << std::endl;
        std::cout << "✓ Metoda Eulera - prostota implementacji" << std::endl;
        std::cout << "✓ Metoda Heuna - ulepszona stabilność" << std::endl;
        std::cout << "✓ Metoda punktu środkowego - dobry kompromis" << std::endl;
        std::cout << "✓ Metoda Runge-Kutta 4 - wysoka dokładność" << std::endl;
        std::cout << "✓ Adaptacyjna RK4/5 - automatyczna kontrola błędu" << std::endl;
        std::cout << "✓ Układy równań różniczkowych" << std::endl;
        std::cout << "✓ Symulacje fizyczne (oscylator, problem trzech ciał)" << std::endl;
        std::cout << "✓ Analiza zbieżności i stabilności" << std::endl;
        
        std::cout << "\nWygenerowane pliki:" << std::endl;
        std::cout << "- ode_basic_example.csv (podstawowe metody)" << std::endl;
        std::cout << "- ode_adaptive_example.csv (metoda adaptacyjna)" << std::endl;
        std::cout << "- oscillator_simulation.csv (symulacja oscylatora)" << std::endl;
        std::cout << "- three_body_simulation.csv (problem trzech ciał)" << std::endl;
        std::cout << "- ode_convergence_analysis.csv (analiza zbieżności)" << std::endl;
        
        std::cout << "\nZastosowania praktyczne:" << std::endl;
        std::cout << "- Dynamika układów mechanicznych" << std::endl;
        std::cout << "- Modelowanie wzrostu populacji" << std::endl;
        std::cout << "- Obwody elektryczne RLC" << std::endl;
        std::cout << "- Reakcje chemiczne" << std::endl;
        std::cout << "- Równania płynów i termodynamiki" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Błąd podczas wykonywania demonstracji: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}