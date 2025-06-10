/**
 * @file nonlinear_equations.h
 * @brief Implementacja metod rozwiązywania równań nieliniowych
 * @details Zawiera algorytmy: bisekcja, Newton-Raphson, sieczne, punkt stały
 */

#ifndef NONLINEAR_EQUATIONS_H
#define NONLINEAR_EQUATIONS_H

#include <vector>
#include <functional>
#include <string>

namespace agh_numerical {

/**
 * @brief Klasa implementująca metody rozwiązywania równań nieliniowych
 */
class NonlinearEquations {
public:
    /**
     * @brief Struktura przechowująca wynik rozwiązywania równania nieliniowego
     */
    struct NonlinearResult {
        double root;                        ///< Znaleziony pierwiastek
        std::vector<double> iterations;     ///< Historia przybliżeń
        std::vector<double> errors;         ///< Historia błędów
        std::vector<double> function_values; ///< Historia wartości funkcji
        int iteration_count;                ///< Liczba iteracji
        bool converged;                     ///< Czy metoda zbiegła
        double final_error;                 ///< Końcowy błąd
        double convergence_rate;            ///< Tempo zbieżności
        std::string method;                 ///< Nazwa użytej metody
        std::string convergence_type;       ///< Typ zbieżności (liniowa/kwadratowa)
        
        NonlinearResult() : root(0.0), iteration_count(0), converged(false), 
                          final_error(0.0), convergence_rate(0.0) {}
    };

    /**
     * @brief Metoda bisekcji dla znajdowania pierwiastków
     * @param f Funkcja, której pierwiastek szukamy
     * @param a Lewa granica przedziału
     * @param b Prawa granica przedziału
     * @param tolerance Tolerancja błędu
     * @param max_iterations Maksymalna liczba iteracji
     * @return Struktura NonlinearResult z wynikiem
     * 
     * @details Metoda bisekcji dzieli przedział na pół w każdej iteracji.
     *          Wymaga f(a)·f(b) < 0 (zmiana znaku na końcach przedziału).
     *          Zbieżność liniowa z tempem ½.
     *          Zawsze zbieżna dla funkcji ciągłych.
     * 
     * @example
     * ```cpp
     * auto f = [](double x) { return x*x - 2.0; }; // sqrt(2)
     * auto result = NonlinearEquations::bisection_method(f, 1.0, 2.0, 1e-8);
     * ```
     */
    static NonlinearResult bisection_method(std::function<double(double)> f,
                                          double a, double b,
                                          double tolerance = 1e-8,
                                          int max_iterations = 100);

    /**
     * @brief Metoda Newtona-Raphsona
     * @param f Funkcja, której pierwiastek szukamy
     * @param df Pochodna funkcji f
     * @param x0 Punkt startowy
     * @param tolerance Tolerancja błędu
     * @param max_iterations Maksymalna liczba iteracji
     * @return Struktura NonlinearResult z wynikiem
     * 
     * @details Metoda Newtona: x_{n+1} = x_n - f(x_n)/f'(x_n)
     *          Zbieżność kwadratowa (przy dobrym punkcie startowym).
     *          Wymaga obliczania pochodnej.
     *          Może nie zbiegać dla złego x0 lub f'(x) ≈ 0.
     */
    static NonlinearResult newton_method(std::function<double(double)> f,
                                       std::function<double(double)> df,
                                       double x0,
                                       double tolerance = 1e-8,
                                       int max_iterations = 50);

    /**
     * @brief Metoda Newtona z numeryczną pochodną
     * @param f Funkcja, której pierwiastek szukamy
     * @param x0 Punkt startowy
     * @param tolerance Tolerancja błędu
     * @param max_iterations Maksymalna liczba iteracji
     * @param h Krok do obliczania pochodnej numerycznej
     * @return Struktura NonlinearResult z wynikiem
     * 
     * @details Pochodna numeryczna: f'(x) ≈ (f(x+h) - f(x-h))/(2h)
     *          Eliminuje potrzebę analitycznego obliczania pochodnej.
     */
    static NonlinearResult newton_numerical_derivative(
        std::function<double(double)> f,
        double x0,
        double tolerance = 1e-8,
        int max_iterations = 50,
        double h = 1e-6);

    /**
     * @brief Metoda siecznych
     * @param f Funkcja, której pierwiastek szukamy
     * @param x0 Pierwszy punkt startowy
     * @param x1 Drugi punkt startowy
     * @param tolerance Tolerancja błędu
     * @param max_iterations Maksymalna liczba iteracji
     * @return Struktura NonlinearResult z wynikiem
     * 
     * @details Metoda siecznych: x_{n+1} = x_n - f(x_n)·(x_n - x_{n-1})/(f(x_n) - f(x_{n-1}))
     *          Aproksymuje pochodną różnicą skończoną.
     *          Zbieżność superlinearna (~1.618).
     *          Nie wymaga obliczania pochodnej.
     */
    static NonlinearResult secant_method(std::function<double(double)> f,
                                       double x0, double x1,
                                       double tolerance = 1e-8,
                                       int max_iterations = 50);

    /**
     * @brief Metoda falsi (regula falsi)
     * @param f Funkcja, której pierwiastek szukamy
     * @param a Lewa granica przedziału
     * @param b Prawa granica przedziału
     * @param tolerance Tolerancja błędu
     * @param max_iterations Maksymalna liczba iteracji
     * @return Struktura NonlinearResult z wynikiem
     * 
     * @details Kombinacja metody bisekcji i siecznych.
     *          Zawsze utrzymuje pierwiastek w przedziale [a,b].
     *          Zbieżność liniowa, ale często szybsza niż bisekcja.
     */
    static NonlinearResult false_position_method(std::function<double(double)> f,
                                                double a, double b,
                                                double tolerance = 1e-8,
                                                int max_iterations = 100);

    /**
     * @brief Metoda iteracji prostej (punkt stały)
     * @param g Funkcja iteracyjna x = g(x)
     * @param x0 Punkt startowy
     * @param tolerance Tolerancja błędu
     * @param max_iterations Maksymalna liczba iteracji
     * @return Struktura NonlinearResult z wynikiem
     * 
     * @details Szuka punktu stałego: x = g(x)
     *          Zbieżność jeśli |g'(x)| < 1 w otoczeniu pierwiastka.
     *          Iteracja: x_{n+1} = g(x_n)
     */
    static NonlinearResult fixed_point_iteration(std::function<double(double)> g,
                                                double x0,
                                                double tolerance = 1e-8,
                                                int max_iterations = 100);

    /**
     * @brief Metoda Brent'a (kombinacja bisekcji i interpolacji)
     * @param f Funkcja, której pierwiastek szukamy
     * @param a Lewa granica przedziału
     * @param b Prawa granica przedziału
     * @param tolerance Tolerancja błędu
     * @return Struktura NonlinearResult z wynikiem
     * 
     * @details Hybrydowa metoda łącząca stabilność bisekcji
     *          z szybkością interpolacji kwadratowej.
     *          Jedna z najnowocześniejszych metod jednowymiarowych.
     */
    static NonlinearResult brent_method(std::function<double(double)> f,
                                      double a, double b,
                                      double tolerance = 1e-12);

    /**
     * @brief Metoda Newtona dla układów równań nieliniowych
     * @param F Wektor funkcji F(x) = 0
     * @param J Macierz Jacobiego (pochodne cząstkowe)
     * @param x0 Wektor punktu startowego
     * @param tolerance Tolerancja błędu
     * @param max_iterations Maksymalna liczba iteracji
     * @return Wektor rozwiązania i informacje o zbieżności
     * 
     * @details Metoda Newtona dla układów: x_{n+1} = x_n - J^{-1}(x_n)·F(x_n)
     *          Wymaga odwracania macierzy Jacobiego w każdej iteracji.
     */
    static std::pair<std::vector<double>, NonlinearResult> newton_system(
        std::function<std::vector<double>(const std::vector<double>&)> F,
        std::function<std::vector<std::vector<double>>(const std::vector<double>&)> J,
        const std::vector<double>& x0,
        double tolerance = 1e-8,
        int max_iterations = 50);

    /**
     * @brief Automatyczne znajdowanie wszystkich pierwiastków w przedziale
     * @param f Funkcja, której pierwiastki szukamy
     * @param a Lewa granica przedziału
     * @param b Prawa granica przedziału
     * @param num_intervals Liczba podprzedziałów do przeszukania
     * @param tolerance Tolerancja dla poszczególnych pierwiastków
     * @return Wektor znalezionych pierwiastków
     * 
     * @details Dzieli przedział [a,b] na podprzedziały i szuka zmian znaku.
     *          Dla każdej zmiany znaku stosuje metodę bisekcji.
     */
    static std::vector<double> find_all_roots(std::function<double(double)> f,
                                            double a, double b,
                                            int num_intervals = 1000,
                                            double tolerance = 1e-8);

    /**
     * @brief Porównanie metod rozwiązywania równań nieliniowych
     * @param f Funkcja testowa
     * @param df Pochodna funkcji (dla metody Newtona)
     * @param a Lewa granica przedziału
     * @param b Prawa granica przedziału
     * @param exact_root Dokładny pierwiastek (do analizy błędów)
     * @return Wektor wyników dla różnych metod
     */
    static std::vector<NonlinearResult> compare_methods(
        std::function<double(double)> f,
        std::function<double(double)> df,
        double a, double b,
        double exact_root = 0.0);

    /**
     * @brief Analiza zbieżności metod numerycznych
     * @param f Funkcja testowa
     * @param df Pochodna funkcji
     * @param exact_root Dokładny pierwiastek
     * @param start_points Wektor punktów startowych do testowania
     * @param filename Plik do zapisu wyników analizy
     */
    static void convergence_analysis(
        std::function<double(double)> f,
        std::function<double(double)> df,
        double exact_root,
        const std::vector<double>& start_points,
        const std::string& filename);

    /**
     * @brief Test stabilności metod dla różnych funkcji
     * @return Raport z testów stabilności
     */
    static std::string stability_test_report();

    /**
     * @brief Rozwiązywanie równania wielomianowego
     * @param coefficients Współczynniki wielomianu (od najwyższego stopnia)
     * @param initial_guesses Początkowe przybliżenia pierwiastków
     * @return Wektor pierwiastków rzeczywistych
     * 
     * @details Używa metody Newtona dla każdego początkowego przybliżenia.
     *          Dla wielomianu P(x) = a_n*x^n + ... + a_1*x + a_0
     */
    static std::vector<double> solve_polynomial(
        const std::vector<double>& coefficients,
        const std::vector<double>& initial_guesses = {});

    /**
     * @brief Metoda optymalizacji (znajdowanie minimum funkcji)
     * @param f Funkcja do minimalizacji
     * @param df Pochodna funkcji
     * @param a Lewa granica przedziału
     * @param b Prawa granica przedziału
     * @param tolerance Tolerancja błędu
     * @return Punkt minimum i wartość funkcji
     * 
     * @details Znajduje minimum poprzez rozwiązanie f'(x) = 0
     *          Używa metody Newtona na pochodnej.
     */
    static std::pair<double, double> find_minimum(
        std::function<double(double)> f,
        std::function<double(double)> df,
        double a, double b,
        double tolerance = 1e-8);

    /**
     * @brief Złoty podział (Golden Section Search) dla optymalizacji
     * @param f Funkcja do minimalizacji (unimodalna)
     * @param a Lewa granica przedziału
     * @param b Prawa granica przedziału
     * @param tolerance Tolerancja błędu
     * @return Punkt minimum i wartość funkcji
     * 
     * @details Metoda nie wymaga pochodnych.
     *          Działa dla funkcji unimodalnych.
     *          Zbieżność liniowa z tempem złotego podziału.
     */
    static std::pair<double, double> golden_section_search(
        std::function<double(double)> f,
        double a, double b,
        double tolerance = 1e-8);

    /**
     * @brief Zapisanie historii zbieżności do pliku CSV
     * @param result Wynik rozwiązywania równania
     * @param filename Nazwa pliku wyjściowego
     * @param include_function_values Czy dołączyć wartości funkcji
     */
    static void save_convergence_history(const NonlinearResult& result,
                                       const std::string& filename,
                                       bool include_function_values = true);

    /**
     * @brief Estymacja tempa zbieżności
     * @param errors Wektor błędów z kolejnych iteracji
     * @return Para [tempo_zbieżności, rząd_zbieżności]
     * 
     * @details Oblicza tempo zbieżności α i rząd p z wzoru:
     *          |e_{n+1}| ≈ α·|e_n|^p
     */
    static std::pair<double, double> estimate_convergence_rate(
        const std::vector<double>& errors);

private:
    /**
     * @brief Obliczenie pochodnej numerycznej (różnice centralne)
     */
    static double numerical_derivative(std::function<double(double)> f, 
                                     double x, double h = 1e-8);

    /**
     * @brief Sprawdzenie warunków zbieżności
     */
    static bool check_convergence(double current, double previous, 
                                double tolerance, int iteration);

    /**
     * @brief Metoda interpolacji kwadratowej (dla metody Brent'a)
     */
    static double quadratic_interpolation(double a, double b, double c,
                                        double fa, double fb, double fc);

    /**
     * @brief Stała złotego podziału
     */
    static const double GOLDEN_RATIO;
    
    static const double EPSILON;           ///< Tolerancja numeryczna
    static const int DEFAULT_MAX_ITER;     ///< Domyślna maksymalna liczba iteracji
};

} // namespace agh_numerical

#endif // NONLINEAR_EQUATIONS_H