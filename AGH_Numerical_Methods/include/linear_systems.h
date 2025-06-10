/**
 * @file linear_systems.h
 * @brief Implementacja metod rozwiązywania układów równań liniowych
 * @details Zawiera algorytmy eliminacji Gaussa z częściowym pivotingiem
 *          oraz rozkład LU dla rozwiązywania układów Ax = b
 */

#ifndef LINEAR_SYSTEMS_H
#define LINEAR_SYSTEMS_H

#include <vector>
#include <string>

namespace agh_numerical {

/**
 * @brief Klasa do rozwiązywania układów równań liniowych
 */
class LinearSystems {
public:
    /**
     * @brief Struktura przechowująca wynik rozwiązania układu równań
     */
    struct Solution {
        std::vector<double> x;      ///< Wektor rozwiązania
        bool is_valid;              ///< Czy rozwiązanie jest poprawne
        double determinant;         ///< Wyznacznik macierzy (jeśli obliczony)
        int iterations;             ///< Liczba iteracji (dla metod iteracyjnych)
        
        Solution() : is_valid(false), determinant(0.0), iterations(0) {}
    };

    /**
     * @brief Rozwiązanie układu równań metodą eliminacji Gaussa z częściowym pivotingiem
     * @param A Macierz współczynników (zostanie zmodyfikowana)
     * @param b Wektor prawych stron (zostanie zmodyfikowany)
     * @return Struktura Solution zawierająca wynik
     * 
     * @details Metoda realizuje eliminację Gaussa z częściowym wyborem elementu podstawowego.
     *          Złożoność obliczeniowa: O(n³)
     * 
     * @example
     * ```cpp
     * std::vector<std::vector<double>> A = {{2, 1}, {1, 3}};
     * std::vector<double> b = {3, 4};
     * auto result = LinearSystems::gauss_elimination(A, b);
     * if (result.is_valid) {
     *     // Rozwiązanie w result.x
     * }
     * ```
     */
    static Solution gauss_elimination(std::vector<std::vector<double>> A, 
                                    std::vector<double> b);

    /**
     * @brief Rozkład LU macierzy z częściowym pivotingiem
     * @param A Macierz wejściowa
     * @param L Macierz dolnotrójkątna (wyjście)
     * @param U Macierz górnotrójkątna (wyjście)
     * @param P Macierz permutacji (wyjście)
     * @return true jeśli rozkład został wykonany pomyślnie
     * 
     * @details Realizuje rozkład PA = LU, gdzie P jest macierzą permutacji
     */
    static bool lu_decomposition(const std::vector<std::vector<double>>& A,
                               std::vector<std::vector<double>>& L,
                               std::vector<std::vector<double>>& U,
                               std::vector<int>& P);

    /**
     * @brief Rozwiązanie układu równań używając rozkładu LU
     * @param L Macierz dolnotrójkątna
     * @param U Macierz górnotrójkątna
     * @param P Wektor permutacji
     * @param b Wektor prawych stron
     * @return Struktura Solution zawierająca wynik
     */
    static Solution solve_lu(const std::vector<std::vector<double>>& L,
                           const std::vector<std::vector<double>>& U,
                           const std::vector<int>& P,
                           const std::vector<double>& b);

    /**
     * @brief Podstawianie w przód dla układu Ly = b
     * @param L Macierz dolnotrójkątna
     * @param b Wektor prawych stron
     * @return Wektor rozwiązania y
     */
    static std::vector<double> forward_substitution(
        const std::vector<std::vector<double>>& L,
        const std::vector<double>& b);

    /**
     * @brief Podstawianie wsteczne dla układu Ux = y
     * @param U Macierz górnotrójkątna
     * @param y Wektor prawych stron
     * @return Wektor rozwiązania x
     */
    static std::vector<double> backward_substitution(
        const std::vector<std::vector<double>>& U,
        const std::vector<double>& y);

    /**
     * @brief Obliczenie wyznacznika macierzy
     * @param A Macierz kwadratowa
     * @return Wartość wyznacznika
     */
    static double determinant(const std::vector<std::vector<double>>& A);

    /**
     * @brief Sprawdzenie poprawności rozwiązania Ax = b
     * @param A Macierz oryginalna
     * @param x Wektor rozwiązania
     * @param b Wektor prawych stron
     * @return Norma błędu residualnego ||Ax - b||
     */
    static double verify_solution(const std::vector<std::vector<double>>& A,
                                const std::vector<double>& x,
                                const std::vector<double>& b);

    /**
     * @brief Wczytanie układu równań z pliku
     * @param filename Ścieżka do pliku
     * @param A Macierz współczynników (wyjście)
     * @param b Wektor prawych stron (wyjście)
     * @return true jeśli wczytano pomyślnie
     * 
     * @details Format pliku:
     * ```
     * N=3
     * b: 1.0 2.0 3.0
     * A:
     * 1.0 2.0 3.0
     * 4.0 5.0 6.0
     * 7.0 8.0 9.0
     * ```
     */
    static bool load_from_file(const std::string& filename,
                             std::vector<std::vector<double>>& A,
                             std::vector<double>& b);

private:
    static const double EPSILON; ///< Tolerancja numeryczna
};

} // namespace agh_numerical

#endif // LINEAR_SYSTEMS_H