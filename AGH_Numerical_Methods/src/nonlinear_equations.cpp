/**
 * @file nonlinear_equations.cpp
 * @brief Implementacja metod rozwiązywania równań nieliniowych
 */

#include "nonlinear_equations.h"
#include "linear_systems.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <chrono>

namespace agh_numerical {

const double NonlinearEquations::GOLDEN_RATIO = (1.0 + std::sqrt(5.0)) / 2.0;
const double NonlinearEquations::EPSILON = 1e-15;
const int NonlinearEquations::DEFAULT_MAX_ITER = 100;

NonlinearEquations::NonlinearResult NonlinearEquations::bisection_method(
    std::function<double(double)> f, double a, double b,
    double tolerance, int max_iterations) {
    
    NonlinearResult result;
    result.method = "Bisection";
    result.convergence_type = "Linear";
    
    double fa = f(a);
    double fb = f(b);
    
    // Sprawdzenie warunków początkowych
    if (fa * fb > 0) {
        std::cerr << "Błąd: Brak zmiany znaku na końcach przedziału [" << a << ", " << b << "]" << std::endl;
        result.converged = false;
        return result;
    }
    
    if (std::abs(fa) < tolerance) {
        result.root = a;
        result.converged = true;
        result.final_error = std::abs(fa);
        return result;
    }
    
    if (std::abs(fb) < tolerance) {
        result.root = b;
        result.converged = true;
        result.final_error = std::abs(fb);
        return result;
    }
    
    double c;
    for (int iter = 0; iter < max_iterations; ++iter) {
        c = (a + b) / 2.0;
        double fc = f(c);
        
        result.iterations.push_back(c);
        result.function_values.push_back(fc);
        
        double error = std::abs(b - a) / 2.0;
        result.errors.push_back(error);
        
        // Sprawdzenie zbieżności
        if (std::abs(fc) < tolerance || error < tolerance) {
            result.root = c;
            result.converged = true;
            result.iteration_count = iter + 1;
            result.final_error = std::min(std::abs(fc), error);
            break;
        }
        
        // Wybór nowego przedziału
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    
    if (!result.converged) {
        result.root = c;
        result.iteration_count = max_iterations;
        result.final_error = std::abs(b - a) / 2.0;
    }
    
    // Obliczenie tempa zbieżności (teoretycznie 0.5 dla bisekcji)
    if (result.errors.size() > 1) {
        result.convergence_rate = 0.5;
    }
    
    return result;
}

NonlinearEquations::NonlinearResult NonlinearEquations::newton_method(
    std::function<double(double)> f, std::function<double(double)> df,
    double x0, double tolerance, int max_iterations) {
    
    NonlinearResult result;
    result.method = "Newton-Raphson";
    result.convergence_type = "Quadratic";
    
    double x = x0;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        double fx = f(x);
        double dfx = df(x);
        
        result.iterations.push_back(x);
        result.function_values.push_back(fx);
        
        // Sprawdzenie czy pochodna nie jest zbyt mała
        if (std::abs(dfx) < EPSILON) {
            std::cerr << "Ostrzeżenie: Pochodna zbyt bliska zeru w punkcie x = " << x << std::endl;
            result.converged = false;
            break;
        }
        
        // Sprawdzenie zbieżności
        if (std::abs(fx) < tolerance) {
            result.root = x;
            result.converged = true;
            result.iteration_count = iter + 1;
            result.final_error = std::abs(fx);
            break;
        }
        
        // Krok Newtona
        double x_new = x - fx / dfx;
        
        double step_error = std::abs(x_new - x);
        result.errors.push_back(step_error);
        
        // Sprawdzenie zbieżności na podstawie zmiany x
        if (step_error < tolerance) {
            result.root = x_new;
            result.converged = true;
            result.iteration_count = iter + 1;
            result.final_error = step_error;
            break;
        }
        
        // Sprawdzenie stabilności
        if (!std::isfinite(x_new)) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w metodzie Newtona" << std::endl;
            result.converged = false;
            break;
        }
        
        x = x_new;
    }
    
    if (!result.converged && result.iteration_count == 0) {
        result.root = x;
        result.iteration_count = max_iterations;
        result.final_error = std::abs(f(x));
    }
    
    // Estymacja tempa zbieżności
    if (result.errors.size() > 2) {
        auto rates = estimate_convergence_rate(result.errors);
        result.convergence_rate = rates.first;
    }
    
    return result;
}

NonlinearEquations::NonlinearResult NonlinearEquations::newton_numerical_derivative(
    std::function<double(double)> f, double x0,
    double tolerance, int max_iterations, double h) {
    
    // Funkcja pochodnej numerycznej
    auto numerical_df = [f, h](double x) {
        return numerical_derivative(f, x, h);
    };
    
    auto result = newton_method(f, numerical_df, x0, tolerance, max_iterations);
    result.method = "Newton (Numerical Derivative)";
    
    return result;
}

NonlinearEquations::NonlinearResult NonlinearEquations::secant_method(
    std::function<double(double)> f, double x0, double x1,
    double tolerance, int max_iterations) {
    
    NonlinearResult result;
    result.method = "Secant";
    result.convergence_type = "Superlinear";
    
    double x_prev = x0;
    double x_curr = x1;
    double f_prev = f(x_prev);
    double f_curr = f(x_curr);
    
    result.iterations.push_back(x_prev);
    result.iterations.push_back(x_curr);
    result.function_values.push_back(f_prev);
    result.function_values.push_back(f_curr);
    
    for (int iter = 1; iter < max_iterations; ++iter) {
        // Sprawdzenie czy różnica wartości funkcji nie jest zbyt mała
        if (std::abs(f_curr - f_prev) < EPSILON) {
            std::cerr << "Ostrzeżenie: Różnica wartości funkcji zbyt mała" << std::endl;
            result.converged = false;
            break;
        }
        
        // Krok metody siecznych
        double x_new = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);
        
        if (!std::isfinite(x_new)) {
            std::cerr << "Ostrzeżenie: Niestabilność numeryczna w metodzie siecznych" << std::endl;
            result.converged = false;
            break;
        }
        
        double f_new = f(x_new);
        
        result.iterations.push_back(x_new);
        result.function_values.push_back(f_new);
        
        double step_error = std::abs(x_new - x_curr);
        result.errors.push_back(step_error);
        
        // Sprawdzenie zbieżności
        if (std::abs(f_new) < tolerance || step_error < tolerance) {
            result.root = x_new;
            result.converged = true;
            result.iteration_count = iter + 1;
            result.final_error = std::min(std::abs(f_new), step_error);
            break;
        }
        
        // Przygotowanie do następnej iteracji
        x_prev = x_curr;
        x_curr = x_new;
        f_prev = f_curr;
        f_curr = f_new;
    }
    
    if (!result.converged) {
        result.root = x_curr;
        result.iteration_count = max_iterations;
        result.final_error = std::abs(f_curr);
    }
    
    // Estymacja tempa zbieżności (teoretycznie ~1.618 dla metody siecznych)
    if (result.errors.size() > 2) {
        auto rates = estimate_convergence_rate(result.errors);
        result.convergence_rate = rates.first;
    }
    
    return result;
}

NonlinearEquations::NonlinearResult NonlinearEquations::false_position_method(
    std::function<double(double)> f, double a, double b,
    double tolerance, int max_iterations) {
    
    NonlinearResult result;
    result.method = "False Position";
    result.convergence_type = "Linear";
    
    double fa = f(a);
    double fb = f(b);
    
    // Sprawdzenie warunków początkowych
    if (fa * fb > 0) {
        std::cerr << "Błąd: Brak zmiany znaku na końcach przedziału" << std::endl;
        result.converged = false;
        return result;
    }
    
    double c;
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Interpolacja liniowa (metoda siecznych z utrzymaniem przedziału)
        c = a - fa * (b - a) / (fb - fa);
        double fc = f(c);
        
        result.iterations.push_back(c);
        result.function_values.push_back(fc);
        
        double error = std::abs(fc);
        result.errors.push_back(error);
        
        // Sprawdzenie zbieżności
        if (error < tolerance) {
            result.root = c;
            result.converged = true;
            result.iteration_count = iter + 1;
            result.final_error = error;
            break;
        }
        
        // Wybór nowego przedziału (utrzymanie zmiany znaku)
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    
    if (!result.converged) {
        result.root = c;
        result.iteration_count = max_iterations;
        result.final_error = std::abs(f(c));
    }
    
    return result;
}

NonlinearEquations::NonlinearResult NonlinearEquations::brent_method(
    std::function<double(double)> f, double a, double b, double tolerance) {
    
    NonlinearResult result;
    result.method = "Brent";
    result.convergence_type = "Superlinear";
    
    double fa = f(a);
    double fb = f(b);
    
    if (fa * fb > 0) {
        std::cerr << "Błąd: Brak zmiany znaku na końcach przedziału" << std::endl;
        result.converged = false;
        return result;
    }
    
    // Upewnij się, że |f(a)| >= |f(b)|
    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    
    double c = a, fc = fa;
    bool mflag = true;
    double s, fs;
    double d = 0;
    
    const int max_iterations = 100;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        result.iterations.push_back(s);
        
        if (std::abs(fb) < tolerance || std::abs(b - a) < tolerance) {
            result.root = b;
            result.converged = true;
            result.iteration_count = iter + 1;
            result.final_error = std::abs(fb);
            break;
        }
        
        if (std::abs(fa - fc) > EPSILON && std::abs(fb - fc) > EPSILON) {
            // Interpolacja kwadratowa
            s = quadratic_interpolation(a, b, c, fa, fb, fc);
        } else {
            // Metoda siecznych
            s = b - fb * (b - a) / (fb - fa);
        }
        
        // Sprawdzenie warunków metody Brenta
        bool condition1 = (s < (3*a + b)/4 || s > b);
        bool condition2 = (mflag && std::abs(s - b) >= std::abs(b - c)/2);
        bool condition3 = (!mflag && std::abs(s - b) >= std::abs(c - d)/2);
        bool condition4 = (mflag && std::abs(b - c) < tolerance);
        bool condition5 = (!mflag && std::abs(c - d) < tolerance);
        
        if (condition1 || condition2 || condition3 || condition4 || condition5) {
            s = (a + b) / 2;
            mflag = true;
        } else {
            mflag = false;
        }
        
        fs = f(s);
        result.function_values.push_back(fs);
        
        d = c;
        c = b;
        fc = fb;
        
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    
    if (!result.converged) {
        result.root = b;
        result.iteration_count = max_iterations;
        result.final_error = std::abs(fb);
    }
    
    return result;
}

std::vector<double> NonlinearEquations::find_all_roots(
    std::function<double(double)> f, double a, double b,
    int num_intervals, double tolerance) {
    
    std::vector<double> roots;
    
    if (num_intervals <= 0 || a >= b) {
        return roots;
    }
    
    double h = (b - a) / num_intervals;
    
    for (int i = 0; i < num_intervals; ++i) {
        double x1 = a + i * h;
        double x2 = a + (i + 1) * h;
        
        double f1 = f(x1);
        double f2 = f(x2);
        
        // Sprawdzenie zmiany znaku
        if (f1 * f2 < 0) {
            auto result = bisection_method(f, x1, x2, tolerance, 50);
            
            if (result.converged) {
                // Sprawdzenie czy pierwiastek już nie został znaleziony
                bool duplicate = false;
                for (double root : roots) {
                    if (std::abs(root - result.root) < tolerance * 10) {
                        duplicate = true;
                        break;
                    }
                }
                
                if (!duplicate) {
                    roots.push_back(result.root);
                }
            }
        }
    }
    
    // Sortowanie pierwiastków
    std::sort(roots.begin(), roots.end());
    
    return roots;
}

std::vector<NonlinearEquations::NonlinearResult> NonlinearEquations::compare_methods(
    std::function<double(double)> f, std::function<double(double)> df,
    double a, double b, double exact_root) {
    
    std::vector<NonlinearResult> results;
    
    // Sprawdzenie czy przedział zawiera pierwiastek
    if (f(a) * f(b) > 0) {
        std::cerr << "Ostrzeżenie: Brak zmiany znaku na końcach przedziału" << std::endl;
    }
    
    // Bisekcja
    auto bisection_result = bisection_method(f, a, b, 1e-8, 100);
    if (bisection_result.converged) {
        results.push_back(bisection_result);
    }
    
    // Newton (punkt startowy w środku przedziału)
    double x0 = (a + b) / 2;
    auto newton_result = newton_method(f, df, x0, 1e-8, 50);
    if (newton_result.converged) {
        results.push_back(newton_result);
    }
    
    // Sieczne
    auto secant_result = secant_method(f, a, b, 1e-8, 50);
    if (secant_result.converged) {
        results.push_back(secant_result);
    }
    
    // Regula falsi
    auto false_pos_result = false_position_method(f, a, b, 1e-8, 100);
    if (false_pos_result.converged) {
        results.push_back(false_pos_result);
    }
    
    // Brent
    auto brent_result = brent_method(f, a, b, 1e-12);
    if (brent_result.converged) {
        results.push_back(brent_result);
    }
    
    // Dodanie informacji o błędzie jeśli znamy dokładny pierwiastek
    if (exact_root != 0.0) {
        for (auto& result : results) {
            double true_error = std::abs(result.root - exact_root);
            result.final_error = true_error;
        }
    }
    
    return results;
}

std::pair<double, double> NonlinearEquations::golden_section_search(
    std::function<double(double)> f, double a, double b, double tolerance) {
    
    double phi = GOLDEN_RATIO;
    double resphi = 2 - phi;
    
    double tol = tolerance;
    double x1 = a + resphi * (b - a);
    double x2 = a + (1 - resphi) * (b - a);
    double f1 = f(x1);
    double f2 = f(x2);
    
    while (std::abs(b - a) > tol) {
        if (f2 > f1) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + (1 - resphi) * (b - a);
            f2 = f(x2);
        }
    }
    
    double x_min = (a + b) / 2;
    return {x_min, f(x_min)};
}

void NonlinearEquations::save_convergence_history(
    const NonlinearResult& result, const std::string& filename,
    bool include_function_values) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Błąd: Nie można utworzyć pliku " << filename << std::endl;
        return;
    }
    
    file << "# Historia zbieżności - metoda: " << result.method << std::endl;
    file << "# Liczba iteracji: " << result.iteration_count << std::endl;
    file << "# Zbieżność: " << (result.converged ? "TAK" : "NIE") << std::endl;
    file << "# Końcowy błąd: " << std::scientific << result.final_error << std::endl;
    
    if (include_function_values && !result.function_values.empty()) {
        file << "iteration,x,f(x),error" << std::endl;
        for (size_t i = 0; i < result.iterations.size(); ++i) {
            file << i << "," << std::fixed << std::setprecision(12) << result.iterations[i];
            
            if (i < result.function_values.size()) {
                file << "," << std::scientific << result.function_values[i];
            } else {
                file << ",";
            }
            
            if (i < result.errors.size()) {
                file << "," << result.errors[i];
            } else {
                file << ",";
            }
            file << std::endl;
        }
    } else {
        file << "iteration,x,error" << std::endl;
        for (size_t i = 0; i < result.iterations.size(); ++i) {
            file << i << "," << std::fixed << std::setprecision(12) << result.iterations[i];
            
            if (i < result.errors.size()) {
                file << "," << std::scientific << result.errors[i];
            } else {
                file << ",";
            }
            file << std::endl;
        }
    }
    
    file.close();
    std::cout << "Historia zbieżności zapisana do pliku: " << filename << std::endl;
}

std::pair<double, double> NonlinearEquations::estimate_convergence_rate(
    const std::vector<double>& errors) {
    
    if (errors.size() < 3) {
        return {1.0, 1.0}; // Zbyt mało danych
    }
    
    // Estymacja tempa zbieżności z ostatnich trzech błędów
    size_t n = errors.size();
    double e_n = errors[n-1];
    double e_n1 = errors[n-2];
    double e_n2 = errors[n-3];
    
    if (e_n <= EPSILON || e_n1 <= EPSILON || e_n2 <= EPSILON) {
        return {1.0, 1.0};
    }
    
    // Oszacowanie rzędu zbieżności p z wzoru: e_{n+1} ≈ C * e_n^p
    double p = std::log(e_n / e_n1) / std::log(e_n1 / e_n2);
    
    // Oszacowanie stałej C
    double C = e_n / std::pow(e_n1, p);
    
    // Ograniczenie do rozsądnych wartości
    p = std::max(0.5, std::min(3.0, p));
    
    return {C, p};
}

double NonlinearEquations::numerical_derivative(
    std::function<double(double)> f, double x, double h) {
    
    // Różnica centralna: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

double NonlinearEquations::quadratic_interpolation(
    double a, double b, double c, double fa, double fb, double fc) {
    
    // Interpolacja kwadratowa przez trzy punkty
    double denom = (fa - fb) * (fa - fc) * (fb - fc);
    
    if (std::abs(denom) < EPSILON) {
        return (a + b) / 2; // Fallback do bisekcji
    }
    
    double numer = a*a*(fb - fc) + b*b*(fc - fa) + c*c*(fa - fb);
    
    return numer / (2 * denom);
}

} // namespace agh_numerical