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
include tests/CMakeFiles/test_linear_systems.dir/depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/test_linear_systems.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/test_linear_systems.dir/flags.make

tests/CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.o: tests/CMakeFiles/test_linear_systems.dir/flags.make
tests/CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.o: ../tests/test_linear_systems.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.o"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.o -c /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/tests/test_linear_systems.cpp

tests/CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.i"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/tests/test_linear_systems.cpp > CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.i

tests/CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.s"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/tests/test_linear_systems.cpp -o CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.s

# Object files for target test_linear_systems
test_linear_systems_OBJECTS = \
"CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.o"

# External object files for target test_linear_systems
test_linear_systems_EXTERNAL_OBJECTS =

tests/test_linear_systems: tests/CMakeFiles/test_linear_systems.dir/test_linear_systems.cpp.o
tests/test_linear_systems: tests/CMakeFiles/test_linear_systems.dir/build.make
tests/test_linear_systems: lib/libagh_numerical_methods.a
tests/test_linear_systems: tests/CMakeFiles/test_linear_systems.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_linear_systems"
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_linear_systems.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/test_linear_systems.dir/build: tests/test_linear_systems

.PHONY : tests/CMakeFiles/test_linear_systems.dir/build

tests/CMakeFiles/test_linear_systems.dir/clean:
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/test_linear_systems.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/test_linear_systems.dir/clean

tests/CMakeFiles/test_linear_systems.dir/depend:
	cd /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/tests /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/tests /workspaces/MN_Projekt_zaliczeniowy/AGH_Numerical_Methods/build/tests/CMakeFiles/test_linear_systems.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/test_linear_systems.dir/depend

