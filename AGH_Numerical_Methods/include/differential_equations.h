/**
 * @file differential_equations.h
 * @brief Implementacja metod rozwiązywania równań różniczkowych zwyczajnych
 * @details Zawiera algorytmy: Euler, Heun, Runge-Kutta, metody adaptacyjne
 */

#ifndef DIFFERENTIAL_EQUATIONS_H
#define DIFFERENTIAL_EQUATIONS_H

#include <vector>
#include <functional>
#include <string>

namespace agh_numerical {

/**
 * @brief Klasa implementująca metody rozwiązywania równań różniczkowych zwyczajnych
 */
class DifferentialEquations {
public:
    /**
     * @brief Typ funkcji prawej strony równania y' = f(t, y)
     */
    using ODEFunction = std::function<double(double, double)>;
    
    /**
     * @brief Typ funkcji prawej strony układu równań y' = f(t, y)
     */
    using ODESystemFunction = std::function<std::vector<double>(double, const std::vector<double>&)>;

    /**
     * @brief Struktura przechowująca rozwiązanie równania różniczkowego
     */
    struct ODESolution {
        std::vector<double> t;          ///< Punkty czasowe
        std::vector<double> y;          ///< Wartości rozwiązania
        std::vector<double> error;      ///< Oszacowanie błędu lokalnego
        double global_error;            ///< Oszacowanie błędu globalnego
        int function_evaluations;      ///< Liczba wywołań funkcji
        bool converged;                 ///< Czy rozwiązanie zbiegło
        std::string method;             ///< Użyta metoda
        
        ODESolution() : global_error(0.0), function_evaluations(0), converged(true) {}
    };

    /**
     * @brief Struktura przechowująca rozwiązanie układu równań różniczkowych
     */
    struct ODESystemSolution {
        std::vector<double> t;                          ///< Punkty czasowe
        std::vector<std::vector<double>> y;             ///< Macierz rozwiązań [czas][zmienna]
        std::vector<std::vector<double>> error;         ///< Oszacowania błędów
        int function_evaluations;                       ///< Liczba wywołań funkcji
        bool converged;                                 ///< Czy rozwiązanie zbiegło
        std::string method;                             ///< Użyta metoda
        
        ODESystemSolution() : function_evaluations(0), converged(true) {}
    };

    /**
     * @brief Metoda Eulera dla równania y' = f(t, y)
     * @param f Funkcja prawej strony równania
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy y(t0) = y0
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @return Struktura ODESolution z rozwiązaniem
     * 
     * @details Metoda Eulera: y_{n+1} = y_n + h·f(t_n, y_n)
     *          Błąd lokalny: O(h²), błąd globalny: O(h)
     *          Metoda jawna pierwszego rzędu.
     * 
     * @example
     * ```cpp
     * auto f = [](double t, double y) { return -y + t; };
     * auto solution = DifferentialEquations::euler_method(f, 0.0, 1.0, 2.0, 0.1);
     * ```
     */
    static ODESolution euler_method(ODEFunction f, 
                                   double t0, double y0, 
                                   double t_end, double h);

    /**
     * @brief Metoda Heuna (ulepszona metoda Eulera)
     * @param f Funkcja prawej strony równania
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @return Struktura ODESolution z rozwiązaniem
     * 
     * @details Metoda Heuna (predyktor-korektor):
     *          Predyktor: ỹ_{n+1} = y_n + h·f(t_n, y_n)
     *          Korektor: y_{n+1} = y_n + (h/2)·[f(t_n, y_n) + f(t_{n+1}, ỹ_{n+1})]
     *          Błąd lokalny: O(h³), błąd globalny: O(h²)
     */
    static ODESolution heun_method(ODEFunction f, 
                                  double t0, double y0, 
                                  double t_end, double h);

    /**
     * @brief Metoda punktu środkowego
     * @param f Funkcja prawej strony równania
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @return Struktura ODESolution z rozwiązaniem
     * 
     * @details Metoda punktu środkowego:
     *          k₁ = f(t_n, y_n)
     *          k₂ = f(t_n + h/2, y_n + h·k₁/2)
     *          y_{n+1} = y_n + h·k₂
     */
    static ODESolution midpoint_method(ODEFunction f, 
                                      double t0, double y0, 
                                      double t_end, double h);

    /**
     * @brief Metoda Rungego-Kutty 4. rzędu (RK4)
     * @param f Funkcja prawej strony równania
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @return Struktura ODESolution z rozwiązaniem
     * 
     * @details Klasyczna metoda RK4:
     *          k₁ = f(t_n, y_n)
     *          k₂ = f(t_n + h/2, y_n + h·k₁/2)
     *          k₃ = f(t_n + h/2, y_n + h·k₂/2)
     *          k₄ = f(t_n + h, y_n + h·k₃)
     *          y_{n+1} = y_n + (h/6)·(k₁ + 2k₂ + 2k₃ + k₄)
     *          Błąd lokalny: O(h⁵), błąd globalny: O(h⁴)
     */
    static ODESolution runge_kutta_4(ODEFunction f, 
                                    double t0, double y0, 
                                    double t_end, double h);

    /**
     * @brief Adaptacyjna metoda Rungego-Kutty z kontrolą błędu
     * @param f Funkcja prawej strony równania
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy
     * @param t_end Końcowy punkt czasowy
     * @param tolerance Żądana dokładność
     * @param h_initial Początkowy krok całkowania
     * @return Struktura ODESolution z rozwiązaniem
     * 
     * @details Używa metod RK4 i RK5 do oszacowania błędu lokalnego
     *          i automatycznego dostosowania kroku całkowania.
     *          Metoda Runge-Kutta-Fehlberg (RKF45).
     */
    static ODESolution adaptive_runge_kutta(ODEFunction f, 
                                          double t0, double y0, 
                                          double t_end, 
                                          double tolerance = 1e-6,
                                          double h_initial = 0.01);

    /**
     * @brief Metoda Eulera dla układu równań różniczkowych
     * @param f Funkcja prawej strony układu
     * @param t0 Początkowy punkt czasowy
     * @param y0 Wektor warunków początkowych
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @return Struktura ODESystemSolution z rozwiązaniem
     */
    static ODESystemSolution euler_system(ODESystemFunction f,
                                         double t0, const std::vector<double>& y0,
                                         double t_end, double h);

    /**
     * @brief Metoda Rungego-Kutty 4. rzędu dla układu równań
     * @param f Funkcja prawej strony układu
     * @param t0 Początkowy punkt czasowy
     * @param y0 Wektor warunków początkowych
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @return Struktura ODESystemSolution z rozwiązaniem
     */
    static ODESystemSolution runge_kutta_4_system(ODESystemFunction f,
                                                  double t0, const std::vector<double>& y0,
                                                  double t_end, double h);

    /**
     * @brief Porównanie metod rozwiązywania równań różniczkowych
     * @param f Funkcja prawej strony równania
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @param analytical_solution Rozwiązanie analityczne (opcjonalne)
     * @return Wektor rozwiązań dla różnych metod
     */
    static std::vector<ODESolution> compare_methods(
        ODEFunction f, double t0, double y0, double t_end, double h,
        std::function<double(double)> analytical_solution = nullptr);

    /**
     * @brief Test stabilności metod dla równania sztywnego
     * @param lambda Parametr sztywności (λ < 0)
     * @param t_end Końcowy czas symulacji
     * @param filename Plik do zapisu wyników
     * @return Analiza stabilności różnych metod
     */
    static void stability_analysis(double lambda, double t_end, 
                                 const std::string& filename);

    /**
     * @brief Rozwiązywanie zagadnienia brzegowego metodą strzałów
     * @param f Funkcja prawej strony równania y'' = f(t, y, y')
     * @param t0 Lewy koniec przedziału
     * @param t1 Prawy koniec przedziału
     * @param y0 Warunek brzegowy y(t0) = y0
     * @param y1 Warunek brzegowy y(t1) = y1
     * @param h Krok całkowania
     * @return Rozwiązanie zagadnienia brzegowego
     */
    static ODESolution boundary_value_shooting_method(
        std::function<double(double, double, double)> f,
        double t0, double t1, double y0, double y1, double h);

    /**
     * @brief Rozwiązywanie równania różniczkowego z opóźnieniem
     * @param f Funkcja prawej strony y'(t) = f(t, y(t), y(t-τ))
     * @param delay Opóźnienie τ
     * @param history Funkcja historii dla t ∈ [t0-τ, t0]
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy
     * @param t_end Końcowy punkt czasowy
     * @param h Krok całkowania
     * @return Struktura ODESolution z rozwiązaniem
     */
    static ODESolution delay_differential_equation(
        std::function<double(double, double, double)> f,
        double delay,
        std::function<double(double)> history,
        double t0, double y0, double t_end, double h);

    /**
     * @brief Analiza zbieżności metod numerycznych
     * @param f Funkcja prawej strony równania
     * @param analytical_solution Dokładne rozwiązanie analityczne
     * @param t0 Początkowy punkt czasowy
     * @param y0 Warunek początkowy
     * @param t_end Końcowy punkt czasowy
     * @param filename Plik do zapisu wyników analizy
     */
    static void convergence_analysis(
        ODEFunction f,
        std::function<double(double)> analytical_solution,
        double t0, double y0, double t_end,
        const std::string& filename);

    /**
     * @brief Symulacja układu oscylatora harmonicznego
     * @param omega Częstość kołowa
     * @param gamma Współczynnik tłumienia
     * @param F0 Amplituda siły wymuszającej
     * @param omega_drive Częstość siły wymuszającej
     * @param x0 Początkowe położenie
     * @param v0 Początkowa prędkość
     * @param t_end Czas symulacji
     * @param h Krok czasowy
     * @return Rozwiązanie układu [pozycja, prędkość]
     * 
     * @details Rozwiązuje równanie: ẍ + 2γẋ + ω²x = F₀cos(ωdt)
     */
    static ODESystemSolution harmonic_oscillator_simulation(
        double omega, double gamma, double F0, double omega_drive,
        double x0, double v0, double t_end, double h);

    /**
     * @brief Symulacja problemu trzech ciał (uproszczona wersja)
     * @param masses Masy trzech ciał
     * @param initial_positions Początkowe pozycje [x1,y1,x2,y2,x3,y3]
     * @param initial_velocities Początkowe prędkości [vx1,vy1,vx2,vy2,vx3,vy3]
     * @param t_end Czas symulacji
     * @param h Krok czasowy
     * @return Rozwiązanie układu dynamicznego
     */
    static ODESystemSolution three_body_problem(
        const std::vector<double>& masses,
        const std::vector<double>& initial_positions,
        const std::vector<double>& initial_velocities,
        double t_end, double h);

    /**
     * @brief Zapisanie rozwiązania do pliku CSV
     * @param solution Rozwiązanie równania różniczkowego
     * @param filename Nazwa pliku wyjściowego
     * @param include_error Czy dołączyć kolumnę z błędem
     */
    static void save_solution_to_csv(const ODESolution& solution,
                                   const std::string& filename,
                                   bool include_error = false);

    /**
     * @brief Zapisanie rozwiązania układu do pliku CSV
     * @param solution Rozwiązanie układu równań
     * @param filename Nazwa pliku wyjściowego
     * @param variable_names Nazwy zmiennych (opcjonalne)
     */
    static void save_system_solution_to_csv(
        const ODESystemSolution& solution,
        const std::string& filename,
        const std::vector<std::string>& variable_names = {});

private:
    /**
     * @brief Funkcja pomocnicza do adaptacyjnej kontroli kroku
     */
    static double compute_optimal_step_size(double h_current, double error, 
                                          double tolerance, int order);

    /**
     * @brief Interpolacja rozwiązania dla równań z opóźnieniem
     */
    static double interpolate_solution(const std::vector<double>& t_hist,
                                     const std::vector<double>& y_hist,
                                     double t_target);

    static const double EPSILON;           ///< Tolerancja numeryczna
    static const double MIN_STEP_SIZE;     ///< Minimalny krok całkowania
    static const double MAX_STEP_SIZE;     ///< Maksymalny krok całkowania
};

} // namespace agh_numerical

#endif // DIFFERENTIAL_EQUATIONS_H