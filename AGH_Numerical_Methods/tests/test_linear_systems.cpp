/**
 * @file test_linear_systems.cpp
 * @brief Testy jednostkowe dla modułu układów równań liniowych
 */

#include "linear_systems.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <functional>
using namespace agh_numerical;

/**
 * @brief Test eliminacji Gaussa dla układu 2x2
 */
bool test_gauss_elimination_2x2() {
    std::cout << "Test: Eliminacja Gaussa 2x2... ";
    
    // Układ: 2x + y = 3, x + 3y = 4 
    // Rozwiązanie: x = 1, y = 1
    std::vector<std::vector<double>> A = {{2.0, 1.0}, {1.0, 3.0}};
    std::vector<double> b = {3.0, 4.0};
    
    auto result = LinearSystems::gauss_elimination(A, b);
    
    if (!result.is_valid) {
        std::cout << "BŁĄD - Nie udało się rozwiązać układu" << std::endl;
        return false;
    }
    
    const double tolerance = 1e-10;
    bool success = (std::abs(result.x[0] - 1.0) < tolerance && 
                   std::abs(result.x[1] - 1.0) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Niepoprawne rozwiązanie: x=" << result.x[0] 
                  << ", y=" << result.x[1] << std::endl;
    }
    
    return success;
}

/**
 * @brief Test eliminacji Gaussa dla układu 3x3
 */
bool test_gauss_elimination_3x3() {
    std::cout << "Test: Eliminacja Gaussa 3x3... ";
    
    // Układ z dokładnym rozwiązaniem x = [1, 2, 3]
    std::vector<std::vector<double>> A = {
        {2.0, 1.0, 1.0},
        {1.0, 3.0, 2.0},
        {1.0, 0.0, 0.0}
    };
    std::vector<double> b = {9.0, 13.0, 1.0};
    
    auto result = LinearSystems::gauss_elimination(A, b);
    
    if (!result.is_valid) {
        std::cout << "BŁĄD - Nie udało się rozwiązać układu" << std::endl;
        return false;
    }
    
    const double tolerance = 1e-10;
    bool success = (std::abs(result.x[0] - 1.0) < tolerance && 
                   std::abs(result.x[1] - 2.0) < tolerance &&
                   std::abs(result.x[2] - 3.0) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Niepoprawne rozwiązanie: ";
        for (size_t i = 0; i < result.x.size(); ++i) {
            std::cout << "x" << i << "=" << result.x[i] << " ";
        }
        std::cout << std::endl;
    }
    
    return success;
}

/**
 * @brief Test rozkładu LU dla macierzy 3x3
 */
bool test_lu_decomposition() {
    std::cout << "Test: Rozkład LU... ";
    
    // Macierz testowa
    std::vector<std::vector<double>> A = {
        {4.0, 3.0, 2.0},
        {3.0, 4.0, 1.0},
        {2.0, 1.0, 3.0}
    };
    
    std::vector<std::vector<double>> L, U;
    std::vector<int> P;
    
    bool decomposition_success = LinearSystems::lu_decomposition(A, L, U, P);
    
    if (!decomposition_success) {
        std::cout << "BŁĄD - Rozkład LU nie powiódł się" << std::endl;
        return false;
    }
    
    // Sprawdzenie wymiarów
    if (L.size() != 3 || U.size() != 3 || P.size() != 3) {
        std::cout << "BŁĄD - Niepoprawne wymiary macierzy L, U lub P" << std::endl;
        return false;
    }
    
    // Sprawdzenie czy L jest dolnotrójkątna z jedynkami na diagonali
    const double tolerance = 1e-10;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(L[i][i] - 1.0) > tolerance) {
            std::cout << "BŁĄD - L nie ma jedynek na diagonali" << std::endl;
            return false;
        }
        for (int j = i + 1; j < 3; ++j) {
            if (std::abs(L[i][j]) > tolerance) {
                std::cout << "BŁĄD - L nie jest dolnotrójkątna" << std::endl;
                return false;
            }
        }
    }
    
    // Sprawdzenie czy U jest górnotrójkątna
    for (int i = 1; i < 3; ++i) {
        for (int j = 0; j < i; ++j) {
            if (std::abs(U[i][j]) > tolerance) {
                std::cout << "BŁĄD - U nie jest górnotrójkątna" << std::endl;
                return false;
            }
        }
    }
    
    std::cout << "POWODZENIE" << std::endl;
    return true;
}

/**
 * @brief Test rozwiązywania układu za pomocą rozkładu LU
 */
bool test_lu_solve() {
    std::cout << "Test: Rozwiązywanie LU... ";
    
    std::vector<std::vector<double>> A = {
        {2.0, 1.0, 1.0},
        {4.0, 3.0, 3.0},
        {8.0, 7.0, 9.0}
    };
    std::vector<double> b = {4.0, 10.0, 24.0};
    
    std::vector<std::vector<double>> L, U;
    std::vector<int> P;
    
    if (!LinearSystems::lu_decomposition(A, L, U, P)) {
        std::cout << "BŁĄD - Rozkład LU nie powiódł się" << std::endl;
        return false;
    }
    
    auto result = LinearSystems::solve_lu(L, U, P, b);
    
    if (!result.is_valid) {
        std::cout << "BŁĄD - Nie udało się rozwiązać układu LU" << std::endl;
        return false;
    }
    
    // Weryfikacja rozwiązania przez podstawienie do oryginalnego układu
    double residual_norm = LinearSystems::verify_solution(A, result.x, b);
    
    const double tolerance = 1e-10;
    bool success = (residual_norm < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE (residual = " << std::scientific 
                  << residual_norm << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - Duży błąd residualny: " << residual_norm << std::endl;
    }
    
    return success;
}

/**
 * @brief Test obliczania wyznacznika
 */
bool test_determinant() {
    std::cout << "Test: Obliczanie wyznacznika... ";
    
    // Macierz 2x2 z wyznacznikiem = 5
    std::vector<std::vector<double>> A2x2 = {{3.0, 1.0}, {2.0, 2.0}};
    double det2x2 = LinearSystems::determinant(A2x2);
    
    // Macierz 3x3 z wyznacznikiem = -2
    std::vector<std::vector<double>> A3x3 = {
        {1.0, 2.0, 3.0},
        {0.0, 1.0, 4.0},
        {5.0, 6.0, 0.0}
    };
    double det3x3 = LinearSystems::determinant(A3x3);
    
    const double tolerance = 1e-10;
    bool success2x2 = (std::abs(det2x2 - 4.0) < tolerance); // 3*2 - 1*2 = 4
    bool success3x3 = (std::abs(det3x3 - 1.0) < tolerance); // Obliczone ręcznie
    
    if (success2x2 && success3x3) {
        std::cout << "POWODZENIE (det2x2=" << det2x2 << ", det3x3=" << det3x3 << ")" << std::endl;
        return true;
    } else {
        std::cout << "BŁĄD - det2x2=" << det2x2 << " (oczek. 4), det3x3=" 
                  << det3x3 << " (oczek. 1)" << std::endl;
        return false;
    }
}

/**
 * @brief Test podstawiania wstecznego
 */
bool test_backward_substitution() {
    std::cout << "Test: Podstawianie wsteczne... ";
    
    // Macierz górnotrójkątna Ux = y
    std::vector<std::vector<double>> U = {
        {2.0, 1.0, 1.0},
        {0.0, 3.0, 2.0},
        {0.0, 0.0, 4.0}
    };
    std::vector<double> y = {8.0, 11.0, 12.0};
    
    auto x = LinearSystems::backward_substitution(U, y);
    
    // Rozwiązanie: x = [1, 2, 3]
    const double tolerance = 1e-10;
    bool success = (x.size() == 3 &&
                   std::abs(x[0] - 1.0) < tolerance &&
                   std::abs(x[1] - 2.0) < tolerance &&
                   std::abs(x[2] - 3.0) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Niepoprawne rozwiązanie: ";
        for (size_t i = 0; i < x.size(); ++i) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }
    
    return success;
}

/**
 * @brief Test podstawiania w przód
 */
bool test_forward_substitution() {
    std::cout << "Test: Podstawianie w przód... ";
    
    // Macierz dolnotrójkątna Ly = b
    std::vector<std::vector<double>> L = {
        {1.0, 0.0, 0.0},
        {2.0, 1.0, 0.0},
        {1.0, 3.0, 1.0}
    };
    std::vector<double> b = {3.0, 8.0, 17.0};
    
    auto y = LinearSystems::forward_substitution(L, b);
    
    // Rozwiązanie: y = [3, 2, 6]
    const double tolerance = 1e-10;
    bool success = (y.size() == 3 &&
                   std::abs(y[0] - 3.0) < tolerance &&
                   std::abs(y[1] - 2.0) < tolerance &&
                   std::abs(y[2] - 6.0) < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE" << std::endl;
    } else {
        std::cout << "BŁĄD - Niepoprawne rozwiązanie: ";
        for (size_t i = 0; i < y.size(); ++i) {
            std::cout << y[i] << " ";
        }
        std::cout << std::endl;
    }
    
    return success;
}

/**
 * @brief Test weryfikacji rozwiązania
 */
bool test_solution_verification() {
    std::cout << "Test: Weryfikacja rozwiązania... ";
    
    std::vector<std::vector<double>> A = {{2.0, 1.0}, {1.0, 3.0}};
    std::vector<double> x = {1.0, 1.0}; // Dokładne rozwiązanie
    std::vector<double> b = {3.0, 4.0};
    
    double residual = LinearSystems::verify_solution(A, x, b);
    
    const double tolerance = 1e-14;
    bool success = (residual < tolerance);
    
    if (success) {
        std::cout << "POWODZENIE (residual = " << std::scientific 
                  << residual << ")" << std::endl;
    } else {
        std::cout << "BŁĄD - Zbyt duży błąd residualny: " << residual << std::endl;
    }
    
    return success;
}

/**
 * @brief Test macierzy osobliwej
 */
bool test_singular_matrix() {
    std::cout << "Test: Macierz osobliwa... ";
    
    // Macierz osobliwa (det = 0)
    std::vector<std::vector<double>> A = {
        {1.0, 2.0, 3.0},
        {2.0, 4.0, 6.0},  // Wiersz 2 = 2 * wiersz 1
        {1.0, 1.0, 1.0}
    };
    std::vector<double> b = {1.0, 2.0, 1.0};
    
    auto result = LinearSystems::gauss_elimination(A, b);
    
    // Oczekujemy, że metoda wykryje osobliwość
    bool success = (!result.is_valid);
    
    if (success) {
        std::cout << "POWODZENIE - Wykryto osobliwość" << std::endl;
    } else {
        std::cout << "BŁĄD - Nie wykryto osobliwości macierzy" << std::endl;
    }
    
    return success;
}

/**
 * @brief Uruchomienie wszystkich testów układów równań liniowych
 */
bool run_all_linear_systems_tests() {
    std::cout << "\n=== TESTY UKŁADÓW RÓWNAŃ LINIOWYCH ===" << std::endl;
    
    std::vector<std::function<bool()>> tests = {
        test_gauss_elimination_2x2,
        test_gauss_elimination_3x3,
        test_lu_decomposition,
        test_lu_solve,
        test_determinant,
        test_backward_substitution,
        test_forward_substitution,
        test_solution_verification,
        test_singular_matrix
    };
    
    int passed = 0;
    int total = tests.size();
    
    for (auto& test : tests) {
        if (test()) {
            passed++;
        }
    }
    
    std::cout << "\nWyniki testów układów równań liniowych: " 
              << passed << "/" << total << " (" 
              << std::fixed << std::setprecision(1) 
              << (100.0 * passed / total) << "%)" << std::endl;
    
    return (passed == total);
}
/**
 * @brief Główna funkcja programu testowego
 */
int main() {
    std::cout << "Uruchamianie testów układów równań liniowych..." << std::endl;
    
    bool all_tests_passed = run_all_linear_systems_tests();
    
    if (all_tests_passed) {
        std::cout << "\n✓ Wszystkie testy przeszły pomyślnie!" << std::endl;
        return 0;
    } else {
        std::cout << "\n✗ Niektóre testy nie przeszły!" << std::endl;
        return 1;
    }
}