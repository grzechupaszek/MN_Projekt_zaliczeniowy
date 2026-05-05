# AGH Numerical Methods Library

**Biblioteka metod numerycznych dla inżynierii obliczeniowej**


## Opis projektu

Biblioteka została zaprojektowana zgodnie z najlepszymi praktykami inżynierii oprogramowania, zapewniając wysoką jakość kodu, dokumentację oraz kompleksowe testowanie funkcjonalności.

## Zakres funkcjonalny

### 🔢 Rozwiązywanie układów równań liniowych
- Eliminacja Gaussa z częściowym pivotingiem
- Rozkład LU z permutacją wierszy
- Obliczanie wyznacznika macierzy
- Weryfikacja poprawności rozwiązań

### 📈 Interpolacja wielomianowa
- Metoda Lagrange'a
- Metoda Newtona z różnicami dzielonymi
- Algorytm Hornera do efektywnego obliczania wielomianów
- Optymalizacja wyboru węzłów interpolacji

### 📊 Aproksymacja funkcji
- Metoda najmniejszych kwadratów
- Aproksymacja wielomianowa
- Analiza błędów aproksymacji

### ∫ Całkowanie numeryczne
- Metoda prostokątów
- Metoda trapezów
- Metoda Simpsona
- Kwadratura Gaussa-Legendre'a
- Metody adaptacyjne

### 🔄 Rozwiązywanie równań różniczkowych
- Metoda Eulera
- Metoda Heuna (Euler ulepszona)
- Metoda Rungego-Kutty 4. rzędu
- Metoda punktu środkowego

### 🎯 Rozwiązywanie równań nieliniowych
- Metoda bisekcji
- Metoda Newtona-Raphsona
- Metoda siecznych
- Analiza zbieżności algorytmów

## Wymagania systemowe

### Minimalne wymagania
- **Kompilator**: GCC 7.0+ lub Clang 6.0+ z obsługą C++17
- **System budowania**: CMake 3.16+
- **System operacyjny**: Linux, macOS, Windows (z MinGW/MSYS2)

### Zalecane środowisko programistyczne
- **IDE**: Visual Studio Code, CLion, Code::Blocks
- **Debugger**: GDB lub LLDB
- **Narzędzia**: Git, Make, Valgrind (Linux)

## Budowanie

### Budowanie biblioteki przy użyciu CMake
```bash
# Tworzenie katalogu dla plików budowania
mkdir build && cd build

# Konfiguracja projektu
cmake .. -DCMAKE_BUILD_TYPE=Release

# Kompilacja
make -j$(nproc)

# Uruchomienie testów jednostkowych
ctest --output-on-failure
```

### Alternatywne budowanie przy użyciu Makefile
```bash
# Kompilacja biblioteki i testów
make all

# Uruchomienie testów
make test

# Czyszczenie plików budowania
make clean
```

### Opcje konfiguracji CMake
- `BUILD_TESTS=ON/OFF` - Włączenie/wyłączenie budowania testów (domyślnie ON)
- `BUILD_EXAMPLES=ON/OFF` - Włączenie/wyłączenie budowania przykładów (domyślnie ON)
- `CMAKE_BUILD_TYPE=Debug/Release` - Tryb budowania (domyślnie Release)

## Struktura projektu

```
AGH_Numerical_Methods/
├── include/                 # Pliki nagłówkowe (.h)
│   ├── numerical_methods.h    # Główny plik nagłówkowy
│   ├── linear_systems.h       # Układy równań liniowych
│   ├── interpolation.h        # Metody interpolacji
│   ├── approximation.h        # Aproksymacja funkcji
│   ├── integration.h          # Całkowanie numeryczne
│   ├── differential_equations.h # Równania różniczkowe
│   └── nonlinear_equations.h   # Równania nieliniowe
├── src/                     # Implementacje (.cpp)
├── tests/                   # Testy jednostkowe
├── examples/               # Przykłady użycia
├── data/                   # Pliki z danymi testowymi
├── docs/                   # Dokumentacja (opcjonalna)
├── CMakeLists.txt         # System budowania CMake
├── Makefile               # Alternatywny system budowania
└── README.md             # Dokumentacja projektu
```

## Przykłady użycia

### Rozwiązywanie układu równań liniowych
```cpp
#include "numerical_methods.h"
using namespace agh_numerical;

int main() {
    // Definicja układu równań: Ax = b
    std::vector<std::vector<double>> A = {
        {2.0, 1.0, -1.0},
        {-3.0, -1.0, 2.0},
        {-2.0, 1.0, 2.0}
    };
    std::vector<double> b = {8.0, -11.0, -3.0};
    
    // Rozwiązanie metodą eliminacji Gaussa
    auto solution = LinearSystems::gauss_elimination(A, b);
    
    if (solution.is_valid) {
        std::cout << "Rozwiązanie: ";
        for (double xi : solution.x) {
            std::cout << xi << " ";
        }
        std::cout << std::endl;
    }
    
    return 0;
}
```

### Interpolacja wielomianowa
```cpp
#include "numerical_methods.h"
using namespace agh_numerical;

int main() {
    // Węzły interpolacji
    std::vector<double> nodes = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> values = {2.0, 4.0, 8.0, 16.0};
    
    // Interpolacja metodą Lagrange'a
    double x = 2.5;
    double result = Interpolation::lagrange_interpolation(nodes, values, x);
    
    std::cout << "f(" << x << ") ≈ " << result << std::endl;
    
    return 0;
}
```

### Całkowanie numeryczne
```cpp
#include "numerical_methods.h"
using namespace agh_numerical;

int main() {
    // Definicja funkcji podcałkowej
    auto f = [](double x) { return x*x + 2*x + 1; };
    
    // Całkowanie metodą Simpsona
    double a = 0.0, b = 2.0;
    int n = 1000;
    double integral = Integration::simpson_method(f, a, b, n);
    
    std::cout << "∫[" << a << "," << b << "] f(x)dx ≈ " << integral << std::endl;
    
    return 0;
}
```

## Testowanie

Biblioteka zawiera kompleksowy zestaw testów jednostkowych dla wszystkich zaimplementowanych algorytmów.

### Uruchomienie testów
```bash
# Za pomocą CMake/CTest
cd build
ctest --verbose

# Bezpośrednie uruchomienie
./run_tests

# Test z dodatkową diagnostyką
make test-verbose
```

### Struktura testów
- **test_linear_systems.cpp** - Testy dla układów równań liniowych
- **test_interpolation.cpp** - Testy metod interpolacji
- **test_integration.cpp** - Testy całkowania numerycznego
- **test_differential_equations.cpp** - Testy równań różniczkowych
- **test_nonlinear_equations.cpp** - Testy równań nieliniowych

## Dokumentacja API

Kompletna dokumentacja API jest dostępna w plikach nagłówkowych w formacie komentarzy Doxygen. Każda funkcja zawiera:
- Opis algorytmu i jego zastosowania
- Specyfikację parametrów wejściowych i wyjściowych
- Przykłady użycia
- Informacje o złożoności obliczeniowej
- Uwagi dotyczące stabilności numerycznej

### Generowanie dokumentacji
```bash
# Wymagane: Doxygen
doxygen Doxyfile
```

## Wydajność i dokładność

### Charakterystyki wydajnościowe
- **Układy równań**: O(n³) dla eliminacji Gaussa, O(n²) dla rozkładu LU
- **Interpolacja**: O(n²) dla Lagrange'a, O(n) dla Hornera
- **Całkowanie**: O(n) dla metod klasycznych, O(1) dla Gaussa-Legendre'a
- **Równania różniczkowe**: O(n) na krok czasowy
- **Równania nieliniowe**: Zbieżność kwadratowa dla metody Newtona

### Dokładność numeryczna
- Wykorzystanie arytmetyki double precision (IEEE 754)
- Implementacja strategii pivotingu dla stabilności numerycznej
- Kontrola błędów zaokrągleń i propagacji błędów
- Walidacja wyników poprzez testy poprawności

## Bibliografia i źródła

1. Burden, R.L., Faires, J.D., Burden, A.M. *Numerical Analysis*, 10th Edition, Cengage Learning, 2015.
2. Press, W.H., Teukolsky, S.A., Vetterling, W.T., Flannery, B.P. *Numerical Recipes in C++*, 3rd Edition, Cambridge University Press, 2007.
3. Kincaid, D., Cheney, W. *Numerical Analysis: Mathematics of Scientific Computing*, 3rd Edition, Brooks/Cole, 2002.
4. Materiały dydaktyczne z przedmiotu "Metody Numeryczne", Wydział Informatyki, Elektroniki i Telekomunikacji AGH.
