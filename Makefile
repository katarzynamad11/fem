# Kompilator
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++17

# Pliki Ÿród³owe
SRCS = Main.cpp Solver.cpp SolveMatrices.cpp GaussQuadrature.cpp DataLoader.cpp Utils.cpp Grid.cpp

# Pliki obiektowe
OBJS = $(SRCS:.cpp=.o)

# Nazwa wyjœciowego programu
TARGET = MES.exe

# Domyœlna regu³a
all: $(TARGET)

# Linkowanie
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Kompilacja pojedynczych plików
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean dla Windows (CMD/PowerShell)
clean:
	del /Q *.o $(TARGET)
