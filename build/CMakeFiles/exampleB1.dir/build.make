# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/abhishek/Project/3_Muoscope

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/abhishek/Project/3_Muoscope/build

# Include any dependencies generated for this target.
include CMakeFiles/exampleB1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/exampleB1.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/exampleB1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/exampleB1.dir/flags.make

CMakeFiles/exampleB1.dir/exampleB1.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/exampleB1.cc.o: ../exampleB1.cc
CMakeFiles/exampleB1.dir/exampleB1.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/exampleB1.dir/exampleB1.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/exampleB1.cc.o -MF CMakeFiles/exampleB1.dir/exampleB1.cc.o.d -o CMakeFiles/exampleB1.dir/exampleB1.cc.o -c /home/abhishek/Project/3_Muoscope/exampleB1.cc

CMakeFiles/exampleB1.dir/exampleB1.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/exampleB1.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/exampleB1.cc > CMakeFiles/exampleB1.dir/exampleB1.cc.i

CMakeFiles/exampleB1.dir/exampleB1.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/exampleB1.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/exampleB1.cc -o CMakeFiles/exampleB1.dir/exampleB1.cc.s

CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o: ../src/ActionInitialization.cc
CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o -MF CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o.d -o CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o -c /home/abhishek/Project/3_Muoscope/src/ActionInitialization.cc

CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/ActionInitialization.cc > CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.i

CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/ActionInitialization.cc -o CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.s

CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o: ../src/DetectorConstruction.cc
CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o -MF CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o.d -o CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o -c /home/abhishek/Project/3_Muoscope/src/DetectorConstruction.cc

CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/DetectorConstruction.cc > CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.i

CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/DetectorConstruction.cc -o CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.s

CMakeFiles/exampleB1.dir/src/EventAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/EventAction.cc.o: ../src/EventAction.cc
CMakeFiles/exampleB1.dir/src/EventAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/exampleB1.dir/src/EventAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/EventAction.cc.o -MF CMakeFiles/exampleB1.dir/src/EventAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/EventAction.cc.o -c /home/abhishek/Project/3_Muoscope/src/EventAction.cc

CMakeFiles/exampleB1.dir/src/EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/EventAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/EventAction.cc > CMakeFiles/exampleB1.dir/src/EventAction.cc.i

CMakeFiles/exampleB1.dir/src/EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/EventAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/EventAction.cc -o CMakeFiles/exampleB1.dir/src/EventAction.cc.s

CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o: ../src/PrimaryGeneratorAction.cc
CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o -MF CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o -c /home/abhishek/Project/3_Muoscope/src/PrimaryGeneratorAction.cc

CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/PrimaryGeneratorAction.cc > CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.i

CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/PrimaryGeneratorAction.cc -o CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.s

CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o: ../src/ROOTManager.cc
CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o -MF CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o.d -o CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o -c /home/abhishek/Project/3_Muoscope/src/ROOTManager.cc

CMakeFiles/exampleB1.dir/src/ROOTManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/ROOTManager.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/ROOTManager.cc > CMakeFiles/exampleB1.dir/src/ROOTManager.cc.i

CMakeFiles/exampleB1.dir/src/ROOTManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/ROOTManager.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/ROOTManager.cc -o CMakeFiles/exampleB1.dir/src/ROOTManager.cc.s

CMakeFiles/exampleB1.dir/src/RunAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/RunAction.cc.o: ../src/RunAction.cc
CMakeFiles/exampleB1.dir/src/RunAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/exampleB1.dir/src/RunAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/RunAction.cc.o -MF CMakeFiles/exampleB1.dir/src/RunAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/RunAction.cc.o -c /home/abhishek/Project/3_Muoscope/src/RunAction.cc

CMakeFiles/exampleB1.dir/src/RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/RunAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/RunAction.cc > CMakeFiles/exampleB1.dir/src/RunAction.cc.i

CMakeFiles/exampleB1.dir/src/RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/RunAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/RunAction.cc -o CMakeFiles/exampleB1.dir/src/RunAction.cc.s

CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o: ../src/SensitiveDetector.cc
CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o -MF CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o.d -o CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o -c /home/abhishek/Project/3_Muoscope/src/SensitiveDetector.cc

CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/SensitiveDetector.cc > CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.i

CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/SensitiveDetector.cc -o CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.s

CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o: ../src/SensitiveDetectorHit.cc
CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o -MF CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o.d -o CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o -c /home/abhishek/Project/3_Muoscope/src/SensitiveDetectorHit.cc

CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/SensitiveDetectorHit.cc > CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.i

CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/SensitiveDetectorHit.cc -o CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.s

CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o: ../src/SteppingAction.cc
CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o -MF CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o -c /home/abhishek/Project/3_Muoscope/src/SteppingAction.cc

CMakeFiles/exampleB1.dir/src/SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/SteppingAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/abhishek/Project/3_Muoscope/src/SteppingAction.cc > CMakeFiles/exampleB1.dir/src/SteppingAction.cc.i

CMakeFiles/exampleB1.dir/src/SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/SteppingAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/abhishek/Project/3_Muoscope/src/SteppingAction.cc -o CMakeFiles/exampleB1.dir/src/SteppingAction.cc.s

# Object files for target exampleB1
exampleB1_OBJECTS = \
"CMakeFiles/exampleB1.dir/exampleB1.cc.o" \
"CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/exampleB1.dir/src/EventAction.cc.o" \
"CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o" \
"CMakeFiles/exampleB1.dir/src/RunAction.cc.o" \
"CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o" \
"CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o" \
"CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o"

# External object files for target exampleB1
exampleB1_EXTERNAL_OBJECTS =

exampleB1: CMakeFiles/exampleB1.dir/exampleB1.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/EventAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/ROOTManager.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/RunAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/SensitiveDetector.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/SensitiveDetectorHit.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/build.make
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4Tree.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4FR.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4GMocren.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4visHepRep.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4RayTracer.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4VRML.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4ToolsSG.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4OpenGL.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4vis_management.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4modeling.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4interfaces.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4persistency.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4error_propagation.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4readout.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4physicslists.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4run.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4event.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4tracking.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4parmodels.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4processes.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4digits_hits.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4track.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4particles.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4geometry.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4materials.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4graphics_reps.so
exampleB1: /home/abhishek/Root/build_root/lib/libCore.so
exampleB1: /home/abhishek/Root/build_root/lib/libImt.so
exampleB1: /home/abhishek/Root/build_root/lib/libRIO.so
exampleB1: /home/abhishek/Root/build_root/lib/libNet.so
exampleB1: /home/abhishek/Root/build_root/lib/libHist.so
exampleB1: /home/abhishek/Root/build_root/lib/libGraf.so
exampleB1: /home/abhishek/Root/build_root/lib/libGraf3d.so
exampleB1: /home/abhishek/Root/build_root/lib/libGpad.so
exampleB1: /home/abhishek/Root/build_root/lib/libROOTDataFrame.so
exampleB1: /home/abhishek/Root/build_root/lib/libTree.so
exampleB1: /home/abhishek/Root/build_root/lib/libTreePlayer.so
exampleB1: /home/abhishek/Root/build_root/lib/libRint.so
exampleB1: /home/abhishek/Root/build_root/lib/libPostscript.so
exampleB1: /home/abhishek/Root/build_root/lib/libMatrix.so
exampleB1: /home/abhishek/Root/build_root/lib/libPhysics.so
exampleB1: /home/abhishek/Root/build_root/lib/libMathCore.so
exampleB1: /home/abhishek/Root/build_root/lib/libThread.so
exampleB1: /home/abhishek/Root/build_root/lib/libMultiProc.so
exampleB1: /home/abhishek/Root/build_root/lib/libROOTVecOps.so
exampleB1: /usr/lib/x86_64-linux-gnu/libGL.so
exampleB1: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.15.3
exampleB1: /usr/lib/x86_64-linux-gnu/libQt5PrintSupport.so.5.15.3
exampleB1: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.15.3
exampleB1: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.15.3
exampleB1: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.15.3
exampleB1: /usr/lib/x86_64-linux-gnu/libXmu.so
exampleB1: /usr/lib/x86_64-linux-gnu/libXext.so
exampleB1: /usr/lib/x86_64-linux-gnu/libXt.so
exampleB1: /usr/lib/x86_64-linux-gnu/libICE.so
exampleB1: /usr/lib/x86_64-linux-gnu/libSM.so
exampleB1: /usr/lib/x86_64-linux-gnu/libX11.so
exampleB1: /usr/lib/x86_64-linux-gnu/libXm.so
exampleB1: /usr/lib/x86_64-linux-gnu/libxerces-c.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4analysis.so
exampleB1: /usr/lib/x86_64-linux-gnu/libexpat.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4zlib.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4intercoms.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4global.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4clhep.so
exampleB1: /home/abhishek/Geant4/geant4_install/lib/libG4ptl.so.2.3.3
exampleB1: CMakeFiles/exampleB1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/abhishek/Project/3_Muoscope/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking CXX executable exampleB1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/exampleB1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/exampleB1.dir/build: exampleB1
.PHONY : CMakeFiles/exampleB1.dir/build

CMakeFiles/exampleB1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/exampleB1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/exampleB1.dir/clean

CMakeFiles/exampleB1.dir/depend:
	cd /home/abhishek/Project/3_Muoscope/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/abhishek/Project/3_Muoscope /home/abhishek/Project/3_Muoscope /home/abhishek/Project/3_Muoscope/build /home/abhishek/Project/3_Muoscope/build /home/abhishek/Project/3_Muoscope/build/CMakeFiles/exampleB1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/exampleB1.dir/depend

