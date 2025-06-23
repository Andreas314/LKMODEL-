# Compiler and flags
CXX = g++
CXXFLAGS = --std=c++23 -O2
LDFLAGS = -llapack -lblas -fopenmp

# Source and object files
SRC = main.cpp \
      BandStructure.cpp \
      Hamiltonian_functions_real.cpp \
      KMesh.cpp \
      Optical_Matrix.cpp \
      Sum.cpp \
      Wrappers.cpp

OBJ = $(SRC:.cpp=.o)

# Executable name
TARGET = lkmodel

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Individual dependencies
main.o: main.cpp Sum.h BandStructure.h KMesh.h Hamiltonian_functions_real.h Wrappers.h Optical_Matrix.h

Sum.o: Sum.cpp Sum.h BandStructure.h KMesh.h Hamiltonian_functions_real.h Wrappers.h Optical_Matrix.h
	$(CXX) $(CXXFLAGS) -fopenmp -c $< -o $@

BandStructure.o: BandStructure.cpp BandStructure.h Hamiltonian_functions_real.h

KMesh.o: KMesh.cpp KMesh.h Hamiltonian_functions_real.h

Hamiltonian_functions_real.o: Hamiltonian_functions_real.cpp Hamiltonian_functions_real.h Wrappers.h

Wrappers.o: Wrappers.cpp Wrappers.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Optical_Matrix.o: Optical_Matrix.cpp Optical_Matrix.h Wrappers.h Hamiltonian_functions_real.h

# Compile all other .cpp files generically
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Cleanup
clean:
	rm -f $(OBJ) $(TARGET)

