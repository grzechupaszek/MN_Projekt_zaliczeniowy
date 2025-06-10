# Makefile dla AGH Numerical Methods Library
# Autor: Student Inżynierii Obliczeniowej AGH
# Alternatywny system budowania dla środowisk bez CMake

# Konfiguracja kompilatora i flag
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic -O2
DEBUG_FLAGS := -g -O0 -DDEBUG
RELEASE_FLAGS := -O3 -DNDEBUG
LDFLAGS := -lm

# Katalogi projektu
INCLUDE_DIR := include
SRC_DIR := src
TEST_DIR := tests
EXAMPLES_DIR := examples
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
BIN_DIR := $(BUILD_DIR)/bin

# Pliki źródłowe i objekty
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
LIBRARY := $(BUILD_DIR)/libagh_numerical_methods.a

# Pliki testowe
TEST_SOURCES := $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS := $(TEST_SOURCES:$(TEST_DIR)/%.cpp=$(OBJ_DIR)/test_%.o)
TEST_BINARY := $(BIN_DIR)/run_tests

# Przykłady użycia
EXAMPLE_SOURCES := $(wildcard $(EXAMPLES_DIR)/*.cpp)
EXAMPLE_BINARIES := $(EXAMPLE_SOURCES:$(EXAMPLES_DIR)/%.cpp=$(BIN_DIR)/%)

# Sprawdzenie trybu budowania
ifdef DEBUG
    CXXFLAGS += $(DEBUG_FLAGS)
    BUILD_TYPE := Debug
else
    CXXFLAGS += $(RELEASE_FLAGS)
    BUILD_TYPE := Release
endif

# Cele główne
.PHONY: all clean test examples install uninstall help info

# Domyślny cel - budowanie biblioteki
all: $(LIBRARY)

# Budowanie biblioteki statycznej
$(LIBRARY): $(OBJECTS) | $(BUILD_DIR)
	@echo "  [AR]  Tworzenie biblioteki statycznej: $@"
	@ar rcs $@ $^
	@echo "  [OK]  Biblioteka utworzona pomyślnie"

# Kompilacja plików objektowych biblioteki
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo "  [CXX] Kompilacja: $< -> $@"
	@$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

# Kompilacja testów jednostkowych
$(OBJ_DIR)/test_%.o: $(TEST_DIR)/%.cpp | $(OBJ_DIR)
	@echo "  [CXX] Kompilacja testu: $< -> $@"
	@$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -I$(TEST_DIR) -c $< -o $@

# Tworzenie pliku wykonywalnego testów
$(TEST_BINARY): $(TEST_OBJECTS) $(LIBRARY) | $(BIN_DIR)
	@echo "  [LD]  Linkowanie testów: $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Kompilacja przykładów użycia
$(BIN_DIR)/%: $(EXAMPLES_DIR)/%.cpp $(LIBRARY) | $(BIN_DIR)
	@echo "  [LD]  Kompilacja przykładu: $< -> $@"
	@$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) $< $(LIBRARY) -o $@ $(LDFLAGS)

# Uruchomienie testów jednostkowych
test: $(TEST_BINARY)
	@echo "  [TEST] Uruchamianie testów jednostkowych..."
	@./$(TEST_BINARY)
	@echo "  [OK]   Testy zakończone"

# Budowanie przykładów użycia
examples: $(EXAMPLE_BINARIES)
	@echo "  [OK]  Wszystkie przykłady skompilowane"

# Tworzenie katalogów budowania
$(BUILD_DIR) $(OBJ_DIR) $(BIN_DIR):
	@echo "  [DIR] Tworzenie katalogu: $@"
	@mkdir -p $@

# Czyszczenie plików budowania
clean:
	@echo "  [RM]  Usuwanie plików budowania..."
	@rm -rf $(BUILD_DIR)
	@echo "  [OK]  Czyszczenie zakończone"

# Instalacja biblioteki w systemie (wymaga uprawnień administratora)
install: $(LIBRARY)
	@echo "  [INSTALL] Instalacja biblioteki w systemie..."
	@sudo mkdir -p /usr/local/include/agh_numerical_methods
	@sudo cp $(INCLUDE_DIR)/*.h /usr/local/include/agh_numerical_methods/
	@sudo cp $(LIBRARY) /usr/local/lib/
	@sudo ldconfig
	@echo "  [OK]     Instalacja zakończona pomyślnie"

# Deinstalacja biblioteki z systemu
uninstall:
	@echo "  [UNINSTALL] Usuwanie biblioteki z systemu..."
	@sudo rm -rf /usr/local/include/agh_numerical_methods
	@sudo rm -f /usr/local/lib/libagh_numerical_methods.a
	@sudo ldconfig
	@echo "  [OK]       Deinstalacja zakończona"

# Informacje o konfiguracji budowania
info:
	@echo "========================================"
	@echo "AGH Numerical Methods Library - Build Info"
	@echo "========================================"
	@echo "Kompilator:       $(CXX)"
	@echo "Flagi kompilacji: $(CXXFLAGS)"
	@echo "Tryb budowania:   $(BUILD_TYPE)"
	@echo "Katalog źródeł:   $(SRC_DIR)"
	@echo "Katalog testów:   $(TEST_DIR)"
	@echo "Katalog budowania: $(BUILD_DIR)"
	@echo "Biblioteka:       $(LIBRARY)"
	@echo "========================================"

# Pomoc - wyświetlenie dostępnych celów
help:
	@echo "AGH Numerical Methods Library - Makefile"
	@echo "========================================="
	@echo "Dostępne cele:"
	@echo "  all        - Budowanie biblioteki statycznej (domyślny)"
	@echo "  test       - Kompilacja i uruchomienie testów jednostkowych"
	@echo "  examples   - Kompilacja przykładów użycia"
	@echo "  clean      - Usunięcie plików budowania"
	@echo "  install    - Instalacja biblioteki w systemie"
	@echo "  uninstall  - Deinstalacja biblioteki z systemu"
	@echo "  info       - Wyświetlenie informacji o konfiguracji"
	@echo "  help       - Wyświetlenie tej pomocy"
	@echo ""
	@echo "Opcje konfiguracji:"
	@echo "  DEBUG=1    - Budowanie w trybie debugowania"
	@echo "  CXX=clang++ - Użycie innego kompilatora"
	@echo ""
	@echo "Przykłady użycia:"
	@echo "  make               # Budowanie biblioteki"
	@echo "  make DEBUG=1       # Budowanie w trybie debug"
	@echo "  make test          # Uruchomienie testów"
	@echo "  make CXX=clang++   # Kompilacja z Clang"
	@echo "  make clean all     # Czyszczenie i ponowne budowanie"

# Cel do pełnego rebuildu
rebuild: clean all

# Cel do uruchomienia konkretnego przykładu
run-example-%: $(BIN_DIR)/%
	@echo "  [RUN]  Uruchamianie przykładu: $*"
	@./$<

# Cel do sprawdzenia składni (kompilacja bez linkowania)
syntax-check:
	@echo "  [CHECK] Sprawdzanie składni plików źródłowych..."
	@for file in $(SOURCES); do \
		echo "    Sprawdzanie: $$file"; \
		$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -fsyntax-only $$file; \
	done
	@echo "  [OK]    Składnia poprawna"

# Cel do analizy statycznej kodu (jeśli cppcheck jest zainstalowany)
static-analysis:
	@echo "  [ANALYSIS] Uruchamianie analizy statycznej..."
	@if command -v cppcheck >/dev/null 2>&1; then \
		cppcheck --enable=all --std=c++17 -I$(INCLUDE_DIR) $(SRC_DIR) $(TEST_DIR); \
	else \
		echo "  [SKIP]     cppcheck nie jest zainstalowany"; \
	fi

# Cel do formatowania kodu (jeśli clang-format jest zainstalowany)
format:
	@echo "  [FORMAT] Formatowanie kodu źródłowego..."
	@if command -v clang-format >/dev/null 2>&1; then \
		find $(SRC_DIR) $(INCLUDE_DIR) $(TEST_DIR) $(EXAMPLES_DIR) -name "*.cpp" -o -name "*.h" | \
		xargs clang-format -i -style="{BasedOnStyle: Google, IndentWidth: 4}"; \
		echo "  [OK]     Formatowanie zakończone"; \
	else \
		echo "  [SKIP]   clang-format nie jest zainstalowany"; \
	fi

# Cel do generowania dokumentacji (jeśli doxygen jest zainstalowany)
docs:
	@echo "  [DOCS]  Generowanie dokumentacji..."
	@if command -v doxygen >/dev/null 2>&1; then \
		if [ -f Doxyfile ]; then \
			doxygen Doxyfile; \
			echo "  [OK]    Dokumentacja wygenerowana w docs/html/"; \
		else \
			echo "  [ERROR] Brak pliku Doxyfile"; \
		fi \
	else \
		echo "  [SKIP]  doxygen nie jest zainstalowany"; \
	fi

# Cel do tworzenia pakietu dystrybucyjnego
package: clean
	@echo "  [PACKAGE] Tworzenie pakietu dystrybucyjnego..."
	@tar -czf AGH_Numerical_Methods_v1.0.0.tar.gz \
		--exclude='.git*' \
		--exclude='*.tar.gz' \
		--exclude='build' \
		.
	@echo "  [OK]      Pakiet utworzony: AGH_Numerical_Methods_v1.0.0.tar.gz"

# Cel do sprawdzenia zależności
check-deps:
	@echo "  [CHECK] Sprawdzanie zależności systemowych..."
	@echo -n "    g++: "; command -v g++ >/dev/null && echo "✓" || echo "✗"
	@echo -n "    make: "; command -v make >/dev/null && echo "✓" || echo "✗"
	@echo -n "    ar: "; command -v ar >/dev/null && echo "✓" || echo "✗"
	@echo "  Opcjonalne narzędzia:"
	@echo -n "    cmake: "; command -v cmake >/dev/null && echo "✓" || echo "✗"
	@echo -n "    cppcheck: "; command -v cppcheck >/dev/null && echo "✓" || echo "✗"
	@echo -n "    clang-format: "; command -v clang-format >/dev/null && echo "✓" || echo "✗"
	@echo -n "    doxygen: "; command -v doxygen >/dev/null && echo "✓" || echo "✗"

# Zapobieganie usuwaniu plików pośrednich
.PRECIOUS: $(OBJ_DIR)/%.o $(OBJ_DIR)/test_%.o

# Deklaracja, że niektóre cele nie tworzą plików
.PHONY: rebuild run-example-% syntax-check static-analysis format docs package check-deps