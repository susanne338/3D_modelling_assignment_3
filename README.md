# 3D_modelling_assignment_3
mainSusanne.cpp is Susanne's file, and main.cpp Puti's/Louis'.

all files submitted:
1. main.cpp -> main c++ file to run
2. json.hpp -> json library file (need this file for main.cpp to run)
3. Duplex_A_20110907.ifc -> the original IFC file
4. Duplex_A_20110907_nofurni.obj -> obj file as the result from IfcConvert
5. report.pdf -> our report
6. README.md

how to run the code
1. All of the input file including: (main.cpp, json.hpp)
must be stored in the same c++ project file
2. The input OBJ file must be put in either debug or release folder, depending in what mode the code will run.
3. This script will gives 2 output in 2 different type (CityJSON and OBJ), both output can be found also in the release/debug folder, depends on what mode the code will run.
4. If the file will be ran multiple times, the OBJ filename must be different than all of the previous one because it will not overwrite.
