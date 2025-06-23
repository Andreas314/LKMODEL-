# Compiler and flags
CXX = g++
CXXFLAGS = --std=c++23 -O2
LDFLAGS = -llapack -lblas -fopenmp

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = .

# Sources and objects
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC))

# Executable name
TARGET = $(BIN_DIR)/lkmodel

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Object file rules with special flags where needed
$(OBJ_DIR)/main.o: $(SRC_DIR)/main.cpp $(SRC_DIR)/Sum.h $(SRC_DIR)/BandStructure.h $(SRC_DIR)/KMesh.h $(SRC_DIR)/Hamiltonian_functions_real.h $(SRC_DIR)/Wrappers.h $(SRC_DIR)/Optical_Matrix.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/Sum.o: $(SRC_DIR)/Sum.cpp $(SRC_DIR)/Sum.h $(SRC_DIR)/BandStructure.h $(SRC_DIR)/KMesh.h $(SRC_DIR)/Hamiltonian_functions_real.h $(SRC_DIR)/Wrappers.h $(SRC_DIR)/Optical_Matrix.h
	$(CXX) $(CXXFLAGS) -fopenmp -c $< -o $@

$(OBJ_DIR)/BandStructure.o: $(SRC_DIR)/BandStructure.cpp $(SRC_DIR)/BandStructure.h $(SRC_DIR)/Hamiltonian_functions_real.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/KMesh.o: $(SRC_DIR)/KMesh.cpp $(SRC_DIR)/KMesh.h $(SRC_DIR)/Hamiltonian_functions_real.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/Hamiltonian_functions_real.o: $(SRC_DIR)/Hamiltonian_functions_real.cpp $(SRC_DIR)/Hamiltonian_functions_real.h $(SRC_DIR)/Wrappers.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/Wrappers.o: $(SRC_DIR)/Wrappers.cpp $(SRC_DIR)/Wrappers.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/Optical_Matrix.o: $(SRC_DIR)/Optical_Matrix.cpp $(SRC_DIR)/Optical_Matrix.h $(SRC_DIR)/Wrappers.h $(SRC_DIR)/Hamiltonian_functions_real.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Generic rule for remaining .cpp files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure object directory exists
$(OBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Cleanup
clean:
	rm -rf $(OBJ_DIR) $(TARGET)
