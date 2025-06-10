/**
 * @file integration.h
 * @brief Implementacja metod całkowania numerycznego
 * @details Zawiera algorytmy całkowania: prostokąty, trapezy, Simpson,
 *          kwadratura Gaussa-Legendre'a oraz metody adaptacyjne
 */

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <vector>
#include <functional>
#include <string>

namespace agh_numerical {

/**
 * @brief Klasa implementująca metody całkowania numerycznego
 */
class Integration {
public:
    /**
     * @brief Struktura przechowująca węzły i wagi kwadratury Gaussa-Legendre'a
     */
    struct GaussLegendreNode {
        double point;   ///< Węzeł kwadratury
        double weight;  ///< Waga węzła
        
        GaussLegendreNode(double p, double w) : point(p), weight(w) {}
    };

    /**
     * @brief Struktura przechowująca wynik całkowania
     */
    struct IntegrationResult {
        double value;           ///< Wartość całki
        double error_estimate;  ///< Oszacowanie błędu
        int function_calls;     ///< Liczba wywołań funkcji
        bool converged;         ///< Czy osiągnięto zbieżność
        std::string method;     ///< Użyta metoda
        
        IntegrationResult() : value(0.0), error_estimate(0.0), 
                            function_calls(0), converged(false) {}
    };

    /**
     * @brief Metoda prostokątów (reguła punktu środkowego)
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param n Liczba podziałów przedziału
     * @return Przybliżona wartość całki
     * 
     * @details Metoda prostokątów przybliża całkę sumą pól prostokątów,
     *          gdzie wysokość każdego prostokąta jest równa wartości funkcji
     *          w punkcie środkowym podprzedziału.
     *          Błąd: O(h²) gdzie h = (b-a)/n
     * 
     * @example
     * ```cpp
     * auto f = [](double x) { return x*x; };
     * double integral = Integration::rectangle_method(f, 0.0, 2.0, 1000);
     * ```
     */
    static double rectangle_method(std::function<double(double)> f,
                                 double a, double b, int n);

    /**
     * @brief Metoda trapezów
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param n Liczba podziałów przedziału
     * @return Przybliżona wartość całki
     * 
     * @details Metoda trapezów przybliża całkę sumą pól trapezów.
     *          Błąd: O(h²) gdzie h = (b-a)/n
     */
    static double trapezoid_method(std::function<double(double)> f,
                                 double a, double b, int n);

    /**
     * @brief Metoda Simpsona (reguła 1/3 Simpsona)
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param n Liczba podziałów przedziału (musi być parzysta)
     * @return Przybliżona wartość całki
     * 
     * @details Metoda Simpsona używa interpolacji wielomianami stopnia 2.
     *          Jest dokładna dla wielomianów stopnia ≤ 3.
     *          Błąd: O(h⁴) gdzie h = (b-a)/n
     */
    static double simpson_method(std::function<double(double)> f,
                               double a, double b, int n);

    /**
     * @brief Pobieranie węzłów i wag kwadratury Gaussa-Legendre'a
     * @param n Liczba węzłów (1-5)
     * @return Wektor węzłów i wag dla przedziału [-1,1]
     * 
     * @details Zwraca predefiniowane węzły i wagi dla kwadratury G-L
     *          na przedziale [-1,1]. Kwadratura n-punktowa jest dokładna
     *          dla wielomianów stopnia ≤ 2n-1.
     */
    static std::vector<GaussLegendreNode> get_gauss_legendre_nodes(int n);

    /**
     * @brief Kwadratura Gaussa-Legendre'a
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param n Liczba węzłów kwadratury (1-5)
     * @return Wartość całki
     * 
     * @details Kwadratura Gaussa-Legendre'a zapewnia najwyższą dokładność
     *          dla danej liczby węzłów. N-punktowa kwadratura jest dokładna
     *          dla wielomianów stopnia ≤ 2N-1.
     */
    static double gauss_legendre_quadrature(std::function<double(double)> f,
                                          double a, double b, int n);

    /**
     * @brief Adaptacyjna kwadratura Gaussa-Legendre'a
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param n Liczba węzłów na podprzedziale
     * @param num_segments Liczba segmentów podziału
     * @return Wartość całki
     * 
     * @details Dzieli przedział całkowania na mniejsze segmenty i stosuje
     *          kwadraturę G-L na każdym z nich. Szczególnie użyteczna
     *          dla funkcji oscylacyjnych lub o szybkim wzroście.
     */
    static double adaptive_gauss_legendre(std::function<double(double)> f,
                                        double a, double b, int n, int num_segments = 10);

    /**
     * @brief Adaptacyjna metoda Simpsona z kontrolą błędu
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param tolerance Żądana dokładność
     * @param max_depth Maksymalna głębokość rekursji
     * @return Struktura IntegrationResult z wynikiem i diagnostyką
     * 
     * @details Rekurencyjnie dzieli przedział na pół, dopóki błąd
     *          nie spadnie poniżej zadanej tolerancji.
     */
    static IntegrationResult adaptive_simpson(std::function<double(double)> f,
                                            double a, double b,
                                            double tolerance = 1e-8,
                                            int max_depth = 15);

    /**
     * @brief Porównanie metod całkowania dla danej funkcji
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param n Liczba podziałów dla metod klasycznych
     * @param exact_value Dokładna wartość całki (dla analizy błędów)
     * @return Wektor wyników dla różnych metod
     */
    static std::vector<IntegrationResult> compare_methods(
        std::function<double(double)> f,
        double a, double b, int n,
        double exact_value = 0.0);

    /**
     * @brief Test zbieżności metody całkowania
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param method Nazwa metody ("rectangle", "trapezoid", "simpson")
     * @param exact_value Dokładna wartość całki
     * @param filename Plik do zapisu wyników
     * @param max_n Maksymalna liczba podziałów
     */
    static void test_convergence(std::function<double(double)> f,
                               double a, double b,
                               const std::string& method,
                               double exact_value,
                               const std::string& filename,
                               int max_n = 4096);

    /**
     * @brief Całkowanie funkcji oscylacyjnych z metodą adaptacyjną
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param estimated_period Szacowany okres oscylacji
     * @param n_per_period Liczba punktów na okres
     * @return Wartość całki
     * 
     * @details Specjalizowana metoda dla funkcji oscylacyjnych,
     *          która automatycznie dostosowuje gęstość siatki do
     *          częstotliwości oscylacji.
     */
    static double integrate_oscillatory(std::function<double(double)> f,
                                      double a, double b,
                                      double estimated_period,
                                      int n_per_period = 10);

    /**
     * @brief Całkowanie funkcji o szybkim wzroście eksponencjalnym
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param growth_parameter Parametr wzrostu (np. dla exp(ax²))
     * @return Wartość całki
     * 
     * @details Używa nietuniformowego podziału z większą gęstością
     *          punktów w obszarach szybkiego wzrostu funkcji.
     */
    static double integrate_exponential_growth(std::function<double(double)> f,
                                             double a, double b,
                                             double growth_parameter);

    /**
     * @brief Pomiar wydajności metod całkowania
     * @param f Funkcja podcałkowa
     * @param a Dolna granica całkowania
     * @param b Górna granica całkowania
     * @param n Liczba podziałów
     * @param iterations Liczba iteracji pomiaru
     * @return Mapa [nazwa_metody -> czas_wykonania_ms]
     */
    static std::vector<std::pair<std::string, double>> benchmark_methods(
        std::function<double(double)> f,
        double a, double b, int n,
        int iterations = 1000);

    /**
     * @brief Sprawdzenie poprawności implementacji na funkcjach testowych
     * @return true jeśli wszystkie testy przeszły pomyślnie
     * 
     * @details Testuje implementację na funkcjach o znanych całkach:
     *          - wielomiany (dokładne dla odpowiednich metod)
     *          - funkcje trygonometryczne
     *          - funkcje wykładnicze
     */
    static bool run_validation_tests();

private:
    /**
     * @brief Funkcja pomocnicza dla adaptacyjnej metody Simpsona
     */
    static IntegrationResult adaptive_simpson_recursive(
        std::function<double(double)> f,
        double a, double b, double tolerance,
        double S, double fa, double fb, double fc,
        int depth, int max_depth, int& function_calls);

    /**
     * @brief Mapowanie punktu z [-1,1] na [a,b]
     */
    static double map_to_interval(double x, double a, double b) {
        return 0.5 * ((b - a) * x + (b + a));
    }

    static const double EPSILON; ///< Tolerancja numeryczna
};

} // namespace agh_numerical

#endif // INTEGRATION_H