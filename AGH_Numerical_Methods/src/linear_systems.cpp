/**
 * @file linear_systems.cpp
 * @brief Implementacja metod rozwiązywania układów równań liniowych
 */

#include "linear_systems.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace agh_numerical {

const double LinearSystems::EPSILON = 1e-12;

LinearSystems::Solution LinearSystems::gauss_elimination(
    std::vector<std::vector<double>> A, 
    std::vector<double> b) {
    
    Solution result;
    int n = A.size();
    
    if (n == 0 || A[0].size() != n || b.size() != n) {
        return result; // Niepoprawne wymiary
    }
    
    // Eliminacja w przód z częściowym pivotingiem
    for (int k = 0; k < n; k++) {
        // Znajdź element główny (pivot)
        int max_row = k;
        double max_val = std::abs(A[k][k]);
        
        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) > max_val) {
                max_val = std::abs(A[i][k]);
                max_row = i;
            }
        }
        
        // Sprawdź czy macierz jest osobliwa
        if (std::abs(A[max_row][k]) < EPSILON) {
            std::cerr << "Ostrzeżenie: Macierz może być osobliwa w kroku " << k << std::endl;
            return result;
        }
        
        // Zamień wiersze jeśli to konieczne
        if (max_row != k) {
            std::swap(A[k], A[max_row]);
            std::swap(b[k], b[max_row]);
        }
        
        // Eliminacja
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    // Podstawianie wsteczne
    result.x = backward_substitution(A, b);
    result.is_valid = true;
    
    return result;
}

bool LinearSystems::lu_decomposition(
    const std::vector<std::vector<double>>& A,
    std::vector<std::vector<double>>& L,
    std::vector<std::vector<double>>& U,
    std::vector<int>& P) {
    
    int n = A.size();
    if (n == 0 || A[0].size() != n) return false;
    
    // Inicjalizacja macierzy
    L.assign(n, std::vector<double>(n, 0.0));
    U = A; // Kopia macierzy A
    P.resize(n);
    
    // Inicjalizacja macierzy permutacji
    for (int i = 0; i < n; i++) {
        P[i] = i;
        L[i][i] = 1.0;
    }
    
    // Rozkład LU z częściowym pivotingiem
    for (int k = 0; k < n; k++) {
        // Znajdź pivot
        int max_row = k;
        double max_val = std::abs(U[k][k]);
        
        for (int i = k + 1; i < n; i++) {
            if (std::abs(U[i][k]) > max_val) {
                max_val = std::abs(U[i][k]);
                max_row = i;
            }
        }
        
        if (std::abs(U[max_row][k]) < EPSILON) {
            return false; // Macierz osobliwa
        }
        
        // Zamień wiersze
        if (max_row != k) {
            std::swap(U[k], U[max_row]);
            std::swap(P[k], P[max_row]);
            
            // Zamień odpowiednie elementy L (tylko dla j < k)
            for (int j = 0; j < k; j++) {
                std::swap(L[k][j], L[max_row][j]);
            }
        }
        
        // Oblicz elementy L i U
        for (int i = k + 1; i < n; i++) {
            L[i][k] = U[i][k] / U[k][k];
            
            for (int j = k; j < n; j++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }
    
    return true;
}

LinearSystems::Solution LinearSystems::solve_lu(
    const std::vector<std::vector<double>>& L,
    const std::vector<std::vector<double>>& U,
    const std::vector<int>& P,
    const std::vector<double>& b) {
    
    Solution result;
    int n = b.size();
    
    // Permutacja wektora b zgodnie z P
    std::vector<double> pb(n);
    for (int i = 0; i < n; i++) {
        pb[i] = b[P[i]];
    }
    
    // Rozwiąż Ly = Pb
    std::vector<double> y = forward_substitution(L, pb);
    
    // Rozwiąż Ux = y
    result.x = backward_substitution(U, y);
    result.is_valid = true;
    
    return result;
}

std::vector<double> LinearSystems::forward_substitution(
    const std::vector<std::vector<double>>& L,
    const std::vector<double>& b) {
    
    int n = b.size();
    std::vector<double> y(n);
    
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
        
        y[i] /= L[i][i];
    }
    
    return y;
}

std::vector<double> LinearSystems::backward_substitution(
    const std::vector<std::vector<double>>& U,
    const std::vector<double>& y) {
    
    int n = y.size();
    std::vector<double> x(n);
    
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        
        x[i] /= U[i][i];
    }
    
    return x;
}

double LinearSystems::determinant(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    if (n == 0 || A[0].size() != n) return 0.0;
    
    std::vector<std::vector<double>> L, U;
    std::vector<int> P;
    
    if (!lu_decomposition(A, L, U, P)) {
        return 0.0; // Macierz osobliwa
    }
    
    // Wyznacznik = iloczyn elementów diagonalnych U * (-1)^(liczba permutacji)
    double det = 1.0;
    int permutation_count = 0;
    
    for (int i = 0; i < n; i++) {
        det *= U[i][i];
        if (P[i] != i) permutation_count++;
    }
    
    return (permutation_count % 2 == 0) ? det : -det;
}

double LinearSystems::verify_solution(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& x,
    const std::vector<double>& b) {
    
    int n = A.size();
    double norm = 0.0;
    
    // Oblicz Ax - b
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        double diff = sum - b[i];
        norm += diff * diff;
    }
    
    return std::sqrt(norm);
}

bool LinearSystems::load_from_file(
    const std::string& filename,
    std::vector<std::vector<double>>& A,
    std::vector<double>& b) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Nie można otworzyć pliku: " << filename << std::endl;
        return false;
    }
    
    std::string line;
    int N = 0;
    
    // Wczytaj rozmiar macierzy
    if (std::getline(file, line)) {
        size_t pos = line.find("N=");
        if (pos != std::string::npos) {
            N = std::stoi(line.substr(pos + 2));
        } else {
            std::cerr << "Błędny format pliku - brak definicji N" << std::endl;
            return false;
        }
    }
    
    if (N <= 0) {
        std::cerr << "Niepoprawny rozmiar macierzy: " << N << std::endl;
        return false;
    }
    
    // Wczytaj wektor b
    if (std::getline(file, line)) {
        if (line.find("b:") != std::string::npos) {
            std::istringstream iss(line.substr(line.find(":") + 1));
            b.clear();
            double val;
            while (iss >> val) {
                b.push_back(val);
            }
            
            if (b.size() != N) {
                std::cerr << "Niepoprawny rozmiar wektora b" << std::endl;
                return false;
            }
        }
    }
    
    // Wczytaj macierz A
    if (std::getline(file, line) && line.find("A:") != std::string::npos) {
        A.assign(N, std::vector<double>(N));
        
        for (int i = 0; i < N; i++) {
            if (std::getline(file, line)) {
                std::istringstream iss(line);
                for (int j = 0; j < N; j++) {
                    if (!(iss >> A[i][j])) {
                        std::cerr << "Błąd odczytu macierzy A[" << i << "][" << j << "]" << std::endl;
                        return false;
                    }
                }
            }
        }
    }
    
    file.close();
    return true;
}

} // namespace agh_numerical