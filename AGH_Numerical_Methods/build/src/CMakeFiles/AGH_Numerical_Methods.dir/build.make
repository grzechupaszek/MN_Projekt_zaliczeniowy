# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build

# Include any dependencies generated for this target.
include src/CMakeFiles/AGH_Numerical_Methods.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/AGH_Numerical_Methods.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make

src/CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.o: src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make
src/CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.o: ../src/approximation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/approximation.cpp

src/CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/approximation.cpp > CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.i

src/CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/approximation.cpp -o CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.s

src/CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.o: src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make
src/CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.o: ../src/differential_equations.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/differential_equations.cpp

src/CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/differential_equations.cpp > CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.i

src/CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/differential_equations.cpp -o CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.s

src/CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.o: src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make
src/CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.o: ../src/integration.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/integration.cpp

src/CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/integration.cpp > CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.i

src/CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/integration.cpp -o CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.s

src/CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.o: src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make
src/CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.o: ../src/interpolation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/interpolation.cpp

src/CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/interpolation.cpp > CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.i

src/CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/interpolation.cpp -o CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.s

src/CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.o: src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make
src/CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.o: ../src/linear_systems.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/linear_systems.cpp

src/CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/linear_systems.cpp > CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.i

src/CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/linear_systems.cpp -o CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.s

src/CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.o: src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make
src/CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.o: ../src/nonlinear_equations.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/nonlinear_equations.cpp

src/CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/nonlinear_equations.cpp > CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.i

src/CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/nonlinear_equations.cpp -o CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.s

src/CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.o: src/CMakeFiles/AGH_Numerical_Methods.dir/flags.make
src/CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.o: ../src/numerical_methods.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/numerical_methods.cpp

src/CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/numerical_methods.cpp > CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.i

src/CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src/numerical_methods.cpp -o CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.s

# Object files for target AGH_Numerical_Methods
AGH_Numerical_Methods_OBJECTS = \
"CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.o" \
"CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.o" \
"CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.o" \
"CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.o" \
"CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.o" \
"CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.o" \
"CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.o"

# External object files for target AGH_Numerical_Methods
AGH_Numerical_Methods_EXTERNAL_OBJECTS =

lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/approximation.cpp.o
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/differential_equations.cpp.o
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/integration.cpp.o
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/interpolation.cpp.o
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/linear_systems.cpp.o
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/nonlinear_equations.cpp.o
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/numerical_methods.cpp.o
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/build.make
lib/libagh_numerical_methods.a: src/CMakeFiles/AGH_Numerical_Methods.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX static library ../lib/libagh_numerical_methods.a"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && $(CMAKE_COMMAND) -P CMakeFiles/AGH_Numerical_Methods.dir/cmake_clean_target.cmake
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/AGH_Numerical_Methods.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/AGH_Numerical_Methods.dir/build: lib/libagh_numerical_methods.a

.PHONY : src/CMakeFiles/AGH_Numerical_Methods.dir/build

src/CMakeFiles/AGH_Numerical_Methods.dir/clean:
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src && $(CMAKE_COMMAND) -P CMakeFiles/AGH_Numerical_Methods.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/AGH_Numerical_Methods.dir/clean

src/CMakeFiles/AGH_Numerical_Methods.dir/depend:
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/src /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/src/CMakeFiles/AGH_Numerical_Methods.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/AGH_Numerical_Methods.dir/depend

