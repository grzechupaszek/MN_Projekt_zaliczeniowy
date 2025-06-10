/**
 * @file approximation.h
 * @brief Implementacja metod aproksymacji funkcji
 * @details Zawiera algorytmy aproksymacji metodą najmniejszych kwadratów,
 *          aproksymację wielomianową, trygonometryczną oraz eksponencjalną
 */

#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#include <vector>
#include <functional>
#include <string>

namespace agh_numerical {

/**
 * @brief Klasa implementująca metody aproksymacji funkcji
 */
class Approximation {
public:
    /**
     * @brief Struktura przechowująca wynik aproksymacji
     */
    struct ApproximationResult {
        std::vector<double> coefficients;  ///< Współczynniki aproksymacji
        double rmse;                       ///< Root Mean Square Error
        double max_error;                  ///< Maksymalny błąd bezwzględny
        double r_squared;                  ///< Współczynnik determinacji R²
        int degree;                        ///< Stopień aproksymacji
        std::string method;                ///< Użyta metoda
        bool is_valid;                     ///< Czy aproksymacja jest poprawna
        
        ApproximationResult() : rmse(0.0), max_error(0.0), r_squared(0.0), 
                              degree(0), is_valid(false) {}
    };

    /**
     * @brief Struktura opisująca funkcję bazową
     */
    struct BasisFunction {
        std::function<double(double)> function;  ///< Funkcja bazowa
        std::string name;                        ///< Nazwa funkcji
        
        BasisFunction(std::function<double(double)> f, const std::string& n) 
            : function(f), name(n) {}
    };

    /**
     * @brief Aproksymacja wielomianowa metodą najmniejszych kwadratów
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x)
     * @param degree Stopień wielomianu aproksymującego
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Znajduje wielomian P(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ
     *          minimalizujący sumę kwadratów błędów Σ(yᵢ - P(xᵢ))²
     *          Używa rozkładu QR lub SVD dla stabilności numerycznej.
     * 
     * @example
     * ```cpp
     * std::vector<double> x = {1, 2, 3, 4, 5};
     * std::vector<double> y = {2.1, 3.9, 6.1, 8.0, 9.9};
     * auto result = Approximation::polynomial_least_squares(x, y, 2);
     * // Aproksymacja wielomianem stopnia 2
     * ```
     */
    static ApproximationResult polynomial_least_squares(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        int degree);

    /**
     * @brief Aproksymacja funkcjami bazowymi (metoda najmniejszych kwadratów)
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x)
     * @param basis_functions Wektor funkcji bazowych
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Ogólna metoda aproksymacji kombinacją liniową funkcji bazowych:
     *          f(x) ≈ c₁φ₁(x) + c₂φ₂(x) + ... + cₙφₙ(x)
     *          gdzie φᵢ(x) są funkcjami bazowymi
     */
    static ApproximationResult basis_function_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        const std::vector<BasisFunction>& basis_functions);

    /**
     * @brief Aproksymacja trygonometryczna (szereg Fouriera dyskretny)
     * @param x_data Wektor punktów x (równomiernie rozłożonych)
     * @param y_data Wektor wartości funkcji f(x)
     * @param num_harmonics Liczba harmonicznych w aproksymacji
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Aproksymacja w postaci:
     *          f(x) ≈ a₀/2 + Σₖ₌₁ⁿ [aₖcos(kωx) + bₖsin(kωx)]
     *          gdzie ω = 2π/(xₘₐₓ - xₘᵢₙ)
     */
    static ApproximationResult trigonometric_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        int num_harmonics);

    /**
     * @brief Aproksymacja eksponencjalna f(x) ≈ a·e^(bx)
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x) (muszą być dodatnie)
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Wykorzystuje linearyzację: ln(f(x)) = ln(a) + bx
     *          Następnie stosuje regresję liniową na przekształconych danych
     */
    static ApproximationResult exponential_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data);

    /**
     * @brief Aproksymacja potęgowa f(x) ≈ a·x^b
     * @param x_data Wektor punktów x (muszą być dodatnie)
     * @param y_data Wektor wartości funkcji f(x) (muszą być dodatnie)
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Wykorzystuje linearyzację: ln(f(x)) = ln(a) + b·ln(x)
     */
    static ApproximationResult power_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data);

    /**
     * @brief Aproksymacja wielomianami Czebyszewa
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x)
     * @param degree Stopień aproksymacji
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Używa wielomianów Czebyszewa jako funkcji bazowych.
     *          Zapewnia lepszą stabilność numeryczną i równomierny błąd aproksymacji.
     */
    static ApproximationResult chebyshev_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        int degree);

    /**
     * @brief Obliczenie wartości aproksymacji w zadanym punkcie
     * @param result Wynik aproksymacji zawierający współczynniki
     * @param x Punkt, w którym obliczamy wartość
     * @return Wartość aproksymacji w punkcie x
     */
    static double evaluate_approximation(const ApproximationResult& result, double x);

    /**
     * @brief Obliczenie statystyk jakości aproksymacji
     * @param x_data Oryginalne punkty x
     * @param y_data Oryginalne wartości y
     * @param result Wynik aproksymacji
     * @return Zaktualizowana struktura result ze statystykami
     */
    static ApproximationResult compute_approximation_statistics(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        ApproximationResult result);

    /**
     * @brief Automatyczny wybór optymalnego stopnia wielomianu
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x)
     * @param max_degree Maksymalny testowany stopień
     * @param validation_split Procent danych do walidacji (0.0-1.0)
     * @return Optimalny stopień wielomianu
     * 
     * @details Używa walidacji krzyżowej lub podziału train/test
     *          do znalezienia stopnia minimalizującego błąd generalizacji
     */
    static int find_optimal_polynomial_degree(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        int max_degree = 10,
        double validation_split = 0.3);

    /**
     * @brief Porównanie różnych metod aproksymacji
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x)
     * @param max_degree Maksymalny stopień dla aproksymacji wielomianowej
     * @return Wektor wyników dla różnych metod
     */
    static std::vector<ApproximationResult> compare_approximation_methods(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        int max_degree = 5);

    /**
     * @brief Aproksymacja funkcji dwóch zmiennych f(x,y) wielomianami
     * @param x_data Wektor współrzędnych x
     * @param y_data Wektor współrzędnych y  
     * @param z_data Wektor wartości funkcji f(x,y)
     * @param degree_x Stopień wielomianu względem x
     * @param degree_y Stopień wielomianu względem y
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Aproksymacja w postaci:
     *          f(x,y) ≈ Σᵢ₌₀ᵐ Σⱼ₌₀ⁿ aᵢⱼ xᵢ yʲ
     */
    static ApproximationResult bivariate_polynomial_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        const std::vector<double>& z_data,
        int degree_x, int degree_y);

    /**
     * @brief Generowanie punktów aproksymacji dla wizualizacji
     * @param result Wynik aproksymacji
     * @param x_min Minimum przedziału
     * @param x_max Maksimum przedziału
     * @param num_points Liczba punktów do wygenerowania
     * @param filename Nazwa pliku wyjściowego
     */
    static void generate_approximation_plot_data(
        const ApproximationResult& result,
        double x_min, double x_max,
        int num_points,
        const std::string& filename);

    /**
     * @brief Aproksymacja funkcji z wagami (ważona metoda najmniejszych kwadratów)
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x)
     * @param weights Wektor wag dla poszczególnych punktów
     * @param degree Stopień wielomianu aproksymującego
     * @return Struktura ApproximationResult z wynikami
     * 
     * @details Minimalizuje ważoną sumę kwadratów: Σwᵢ(yᵢ - P(xᵢ))²
     */
    static ApproximationResult weighted_polynomial_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        const std::vector<double>& weights,
        int degree);

    /**
     * @brief Aproksymacja robustowa (odporna na wartości odstające)
     * @param x_data Wektor punktów x
     * @param y_data Wektor wartości funkcji f(x)
     * @param degree Stopień wielomianu aproksymującego
     * @param method Metoda robustowa ("huber", "bisquare", "cauchy")
     * @return Struktura ApproximationResult z wynikami
     */
    static ApproximationResult robust_approximation(
        const std::vector<double>& x_data,
        const std::vector<double>& y_data,
        int degree,
        const std::string& method = "huber");

private:
    /**
     * @brief Rozwiązanie układu równań normalnych metodą Choleskiego
     * @param A Macierz A^T*A
     * @param b Wektor A^T*y
     * @return Wektor współczynników
     */
    static std::vector<double> solve_normal_equations(
        const std::vector<std::vector<double>>& A,
        const std::vector<double>& b);

    /**
     * @brief Obliczenie wielomianu Czebyszewa n-tego stopnia
     * @param n Stopień wielomianu
     * @param x Argument
     * @return Wartość Tₙ(x)
     */
    static double chebyshev_polynomial(int n, double x);

    /**
     * @brief Mapowanie z przedziału [a,b] na [-1,1] dla wielomianów Czebyszewa
     */
    static double map_to_chebyshev_interval(double x, double a, double b);

    /**
     * @brief Funkcja wagi dla aproksymacji robustowej
     */
    static double robust_weight_function(double residual, const std::string& method);

    static const double EPSILON; ///< Tolerancja numeryczna
};

} // namespace agh_numerical

#endif // APPROXIMATION_H