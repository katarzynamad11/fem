#pragma once
#include <vector> 


const int npc = 4;
const int num = 2;

struct Node {
    double x, y;
    int BC;  // Boundary condition
};

struct Element {
    int* ID = new int[4];
    double HBC[4][4];
    double P[4]; // wektor 
    double C_local[4][4];
};

struct Grid {
    int nN;     // liczba wezlow - 16 
    int nE;     // liczba elementow - 9 
    std::vector<Node> ND; // wektor przechowuj¹cy wezly
    std::vector<Element> El; // wektor przechowuj¹cy elementy 
};

struct GlobalData {
    double simulationTime;
    double simulationStepTime;
    double conductivity;
    double alfa;
    double tot;
    double initialTemp;
    double density;
    double specificHeat;
};

struct ElemUniv {
    std::vector<double> ksi;      // Wartoœci ksi dla punktów ca³kowania
    std::vector<double> eta;      // Wartoœci eta dla punktów ca³kowania
    std::vector<std::vector<double>> dN_dKsi; // Pochodne dN/dksi
    std::vector<std::vector<double>> dN_dEta; // Pochodne dN/deta
    std::vector<double> w1;      // Wartoœci ksi dla punktów ca³kowania


};

struct Jakobian {
    double J[2][2][npc];
    double detJ[npc];
    double J_inv[2][2][npc];

    std::vector<std::vector<double>> dN_dX;
    std::vector<std::vector<double>> dN_dY;

    Jakobian() {
        dN_dX.resize(npc, std::vector<double>(4, 0.0));
        dN_dY.resize(npc, std::vector<double>(4, 0.0));
    }
};

struct H {
    double Hpcx[4][4][npc];  // Macierz Hpcx dla ka¿dego punktu Gaussa (4x4)
    double Hpcy[4][4][npc];  // Macierz Hpcy dla ka¿dego punktu Gaussa (4x4)
    double Hpcxy[4][4][npc]; // Macierz Hpcxy (suma Hpcx i Hpcy)
    double Hpc[4][4][npc];   // Macierz Hpc (wynik mno¿enia k * detJ * Hpcxy)

    double H_lok[4][4];
    double Hbc[4][4];
    std::vector<std::vector<double>> H_globalna;
    std::vector<std::vector<std::vector<double>>> H_matrix;
    H(int numElements, int globalSize) {
        H_matrix.resize(numElements, std::vector<std::vector<double>>(4, std::vector<double>(4, 0.0)));
        H_globalna.resize(globalSize, std::vector<double>(globalSize, 0.0));
    }
    H() {
        // Inicjalizacja macierzy 4x4 dla ka¿dego punktu Gaussa
        for (int pc = 0; pc < npc; ++pc) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    Hpcx[i][j][pc] = 0.0;
                    Hpcy[i][j][pc] = 0.0;
                    Hpcxy[i][j][pc] = 0.0;
                    Hpc[i][j][pc] = 0.0;
                }
            }
        }
    }
};

struct elementC {
    std::vector<std::vector<double>> N;  // Macierz funkcji kszta³tu 
    std::vector<std::vector<double>> Cpc;  // Macierz Cpc
};