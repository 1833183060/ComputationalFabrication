# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /home/aespielberg/Downloads/CLion-2018.2.2/clion-2018.2.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/aespielberg/Downloads/CLion-2018.2.2/clion-2018.2.2/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aespielberg/ResearchCode/OpenFab

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug

# Include any dependencies generated for this target.
include ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/depend.make

# Include the progress variables for this target.
include ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/progress.make

# Include the compile flags for this target's objects.
include ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/flags.make

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/tearing.c.o: ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/flags.make
ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/tearing.c.o: ../ExternalLibs/libigl/external/glfw/tests/tearing.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/tearing.c.o"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/tearing.dir/tearing.c.o   -c /home/aespielberg/ResearchCode/OpenFab/ExternalLibs/libigl/external/glfw/tests/tearing.c

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/tearing.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/tearing.c.i"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/aespielberg/ResearchCode/OpenFab/ExternalLibs/libigl/external/glfw/tests/tearing.c > CMakeFiles/tearing.dir/tearing.c.i

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/tearing.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/tearing.c.s"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/aespielberg/ResearchCode/OpenFab/ExternalLibs/libigl/external/glfw/tests/tearing.c -o CMakeFiles/tearing.dir/tearing.c.s

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o: ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/flags.make
ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o: ../ExternalLibs/libigl/external/glfw/deps/glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/tearing.dir/__/deps/glad.c.o   -c /home/aespielberg/ResearchCode/OpenFab/ExternalLibs/libigl/external/glfw/deps/glad.c

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/__/deps/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/tearing.dir/__/deps/glad.c.i"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/aespielberg/ResearchCode/OpenFab/ExternalLibs/libigl/external/glfw/deps/glad.c > CMakeFiles/tearing.dir/__/deps/glad.c.i

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/__/deps/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/tearing.dir/__/deps/glad.c.s"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/aespielberg/ResearchCode/OpenFab/ExternalLibs/libigl/external/glfw/deps/glad.c -o CMakeFiles/tearing.dir/__/deps/glad.c.s

# Object files for target tearing
tearing_OBJECTS = \
"CMakeFiles/tearing.dir/tearing.c.o" \
"CMakeFiles/tearing.dir/__/deps/glad.c.o"

# External object files for target tearing
tearing_EXTERNAL_OBJECTS =

ExternalLibs/libigl/external/glfw/tests/tearing: ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/tearing.c.o
ExternalLibs/libigl/external/glfw/tests/tearing: ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/__/deps/glad.c.o
ExternalLibs/libigl/external/glfw/tests/tearing: ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/build.make
ExternalLibs/libigl/external/glfw/tests/tearing: ExternalLibs/libigl/external/glfw/src/libglfw3.a
ExternalLibs/libigl/external/glfw/tests/tearing: /usr/lib/x86_64-linux-gnu/libm.so
ExternalLibs/libigl/external/glfw/tests/tearing: /usr/lib/x86_64-linux-gnu/librt.so
ExternalLibs/libigl/external/glfw/tests/tearing: /usr/lib/x86_64-linux-gnu/libm.so
ExternalLibs/libigl/external/glfw/tests/tearing: /usr/lib/x86_64-linux-gnu/libX11.so
ExternalLibs/libigl/external/glfw/tests/tearing: ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable tearing"
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tearing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/build: ExternalLibs/libigl/external/glfw/tests/tearing

.PHONY : ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/build

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/clean:
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests && $(CMAKE_COMMAND) -P CMakeFiles/tearing.dir/cmake_clean.cmake
.PHONY : ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/clean

ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/depend:
	cd /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aespielberg/ResearchCode/OpenFab /home/aespielberg/ResearchCode/OpenFab/ExternalLibs/libigl/external/glfw/tests /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests /home/aespielberg/ResearchCode/OpenFab/cmake-build-debug/ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ExternalLibs/libigl/external/glfw/tests/CMakeFiles/tearing.dir/depend

