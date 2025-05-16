# Kompilator
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

# Pliki �r�d�owe
SRCS = Main.cpp Solver.cpp SolveMatrices.cpp GaussQuadrature.cpp DataLoader.cpp Utils.cpp Grid.cpp

# Pliki obiektowe
OBJS = $(SRCS:.cpp=.o)

# Nazwa wyj�ciowego programu
TARGET = MES.exe

# Domy�lna regu�a
all: $(TARGET)

# Linkowanie
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Kompilacja pojedynczych plik�w
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean dla Windows (CMD/PowerShell)
clean:
	del /Q *.o $(TARGET)
