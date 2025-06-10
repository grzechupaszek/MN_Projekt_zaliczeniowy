/**
 * @file interpolation.h
 * @brief Implementacja metod interpolacji wielomianowej
 * @details Zawiera algorytmy interpolacji Lagrange'a, Newtona oraz
 *          algorytm Hornera do efektywnego obliczania wartości wielomianów
 */

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <string>
#include <functional>

namespace agh_numerical {

/**
 * @brief Klasa implementująca metody interpolacji wielomianowej
 */
class Interpolation {
public:
    /**
     * @brief Struktura przechowująca wynik interpolacji
     */
    struct InterpolationResult {
        std::vector<double> coefficients;  ///< Współczynniki wielomianu
        std::vector<double> nodes;         ///< Węzły interpolacji
        std::vector<double> values;        ///< Wartości w węzłach
        double mse;                        ///< Średni błąd kwadratowy
        bool is_valid;                     ///< Czy interpolacja jest poprawna
        
        InterpolationResult() : mse(0.0), is_valid(false) {}
    };

    /**
     * @brief Interpolacja wielomianowa metodą Lagrange'a
     * @param nodes Węzły interpolacji (x_i)
     * @param values Wartości funkcji w węzłach (f(x_i))
     * @param x Punkt, w którym obliczamy wartość interpolacji
     * @return Wartość wielomianu interpolacyjnego w punkcie x
     * 
     * @details Metoda Lagrange'a konstruuje wielomian interpolacyjny poprzez
     *          sumę iloczynów wartości funkcji i odpowiadających im
     *          wielomianów bazowych Lagrange'a.
     *          Złożoność obliczeniowa: O(n²)
     * 
     * @example
     * ```cpp
     * std::vector<double> x = {1.0, 2.0, 3.0};
     * std::vector<double> y = {2.0, 4.0, 8.0};
     * double result = Interpolation::lagrange_interpolation(x, y, 2.5);
     * ```
     */
    static double lagrange_interpolation(const std::vector<double>& nodes,
                                       const std::vector<double>& values,
                                       double x);

    /**
     * @brief Obliczenie tablicy różnic dzielonych dla interpolacji Newtona
     * @param nodes Węzły interpolacji
     * @param values Wartości funkcji w węzłach
     * @return Tablica różnic dzielonych
     * 
     * @details Różnice dzielone są podstawą dla konstrukcji wielomianu
     *          interpolacyjnego w postaci Newtona. Tablica jest wypełniana
     *          zgodnie ze wzorem rekurencyjnym.
     */
    static std::vector<std::vector<double>> compute_divided_differences(
        const std::vector<double>& nodes,
        const std::vector<double>& values);

    /**
     * @brief Interpolacja wielomianowa metodą Newtona
     * @param nodes Węzły interpolacji
     * @param divided_diffs Tablica różnic dzielonych
     * @param x Punkt obliczenia
     * @return Wartość wielomianu w punkcie x
     * 
     * @details Wykorzystuje postać Newtona wielomianu interpolacyjnego
     *          z różnicami dzielonymi. Bardziej efektywna od metody
     *          Lagrange'a przy wielokrotnych obliczeniach.
     */
    static double newton_interpolation(const std::vector<double>& nodes,
                                     const std::vector<std::vector<double>>& divided_diffs,
                                     double x);

    /**
     * @brief Algorytm Hornera do obliczania wartości wielomianu
     * @param coefficients Współczynniki wielomianu w porządku malejącym
     * @param x Punkt obliczenia
     * @return Wartość wielomianu P(x)
     * 
     * @details Schemat Hornera pozwala na efektywne obliczenie wartości
     *          wielomianu przy użyciu minimalnej liczby mnożeń.
     *          Dla wielomianu a_n*x^n + ... + a_1*x + a_0
     *          Złożoność: O(n) mnożeń zamiast O(n²)
     * 
     * @example
     * ```cpp
     * // Wielomian 2x³ + 3x² - x + 5
     * std::vector<double> coeffs = {2, 3, -1, 5};
     * double result = Interpolation::horner_evaluation(coeffs, 2.0);
     * ```
     */
    static double horner_evaluation(const std::vector<double>& coefficients,
                                  double x);

    /**
     * @brief Forma naturalna obliczania wielomianu (dla porównania z Hornerem)
     * @param coefficients Współczynniki wielomianu
     * @param x Punkt obliczenia
     * @return Wartość wielomianu
     */
    static double natural_form_evaluation(const std::vector<double>& coefficients,
                                        double x);

    /**
     * @brief Wybór optymalnych węzłów interpolacji co k-ty punkt
     * @param all_nodes Wszystkie dostępne węzły
     * @param all_values Wszystkie wartości funkcji
     * @param step Krok wyboru węzłów (co step-ty punkt)
     * @param selected_nodes Wybrane węzły (wyjście)
     * @param selected_values Wybrane wartości (wyjście)
     */
    static void select_nodes(const std::vector<double>& all_nodes,
                           const std::vector<double>& all_values,
                           int step,
                           std::vector<double>& selected_nodes,
                           std::vector<double>& selected_values);

    /**
     * @brief Obliczenie średniego błędu kwadratowego interpolacji
     * @param all_nodes Wszystkie punkty testowe
     * @param all_values Rzeczywiste wartości funkcji
     * @param interp_nodes Węzły użyte do interpolacji
     * @param interp_values Wartości w węzłach interpolacji
     * @return Wartość MSE dla punktów nie będących węzłami
     * 
     * @details MSE = sqrt(1/n * Σ(y_i - ŷ_i)²) gdzie y_i to rzeczywiste
     *          wartości, a ŷ_i to wartości interpolowane
     */
    static double compute_mse(const std::vector<double>& all_nodes,
                            const std::vector<double>& all_values,
                            const std::vector<double>& interp_nodes,
                            const std::vector<double>& interp_values);

    /**
     * @brief Znajdowanie optymalnej liczby węzłów interpolacji
     * @param nodes Wszystkie dostępne węzły
     * @param values Wszystkie wartości
     * @param min_nodes Minimalna liczba węzłów do testowania
     * @param max_nodes Maksymalna liczba węzłów do testowania
     * @return Optymalna liczba węzłów (minimalizująca MSE)
     */
    static int find_optimal_nodes_count(const std::vector<double>& nodes,
                                      const std::vector<double>& values,
                                      int min_nodes = 3,
                                      int max_nodes = 20);

    /**
     * @brief Wczytanie danych do interpolacji z pliku
     * @param filename Ścieżka do pliku z danymi
     * @param nodes Węzły interpolacji (wyjście)
     * @param values Wartości funkcji (wyjście)
     * @return true jeśli wczytano pomyślnie
     * 
     * @details Format pliku:
     * ```
     * xi: 1.0 2.0 3.0 4.0
     * f(xi): 2.0 4.0 8.0 16.0
     * ```
     */
    static bool load_data(const std::string& filename,
                        std::vector<double>& nodes,
                        std::vector<double>& values);

    /**
     * @brief Generowanie danych interpolacji do pliku (dla wizualizacji)
     * @param nodes Węzły interpolacji
     * @param values Wartości w węzłach
     * @param filename Nazwa pliku wyjściowego
     * @param num_points Liczba punktów do wygenerowania
     * @param method Metoda interpolacji ("lagrange" lub "newton")
     */
    static void generate_interpolation_data(const std::vector<double>& nodes,
                                          const std::vector<double>& values,
                                          const std::string& filename,
                                          int num_points = 1000,
                                          const std::string& method = "lagrange");

    /**
     * @brief Porównanie wydajności metod Hornera i formy naturalnej
     * @param coefficients Współczynniki wielomianu
     * @param x Punkt testowy
     * @param iterations Liczba iteracji testu
     * @return Para (czas_horner, czas_naturalna) w sekundach
     */
    static std::pair<double, double> compare_evaluation_methods(
        const std::vector<double>& coefficients,
        double x,
        int iterations = 1000000);

private:
    static const double EPSILON; ///< Tolerancja numeryczna
};

} // namespace agh_numerical

#endif // INTERPOLATION_H