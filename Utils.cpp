#include <iostream>
#include <iomanip>  
#include <vector>
#include "DataLoader.h"
#include "Grid.h"



void printElements(const Grid& grid) {
    for (int i = 0; i < grid.nE; ++i) {
        std::cout << "Element ID: " << i + 1 << std::endl;
        std::cout << "Nodes: ";
        for (int j = 0; j < 4; ++j) {
            int node_id = grid.El[i].ID[j] - 1; // Indeksowanie wêz³ów od 1 w pliku, wiêc odejmujemy 1
            const Node& node = grid.ND[node_id];
            std::cout << "(" << node.x << ", " << node.y << ") ";
        }
        std::cout << std::endl;
    }
    for (size_t i = 0; i < grid.ND.size(); ++i) {
        const Node& node = grid.ND[i];  // Dostêp do konkretnego wêz³a
        std::cout << "Node ID: " << i + 1 << " BC: " << node.BC << std::endl;
    }
}

void printDerivatives(const ElemUniv& elem, int npc) {


    std::cout << "Integration points (N=" << npc << "):" << std::endl;
    for (size_t i = 0; i < elem.ksi.size(); ++i) {
        std::cout << "Point " << i + 1 << ": ksi = " << elem.ksi[i] << ", eta = " << elem.eta[i] << std::endl;
    }

    std::cout << "\nPochodne funkcji ksztaltu po ksi:" << std::endl;
    std::cout << " pc     dN1/dksi    dN2/dksi    dN3/dksi    dN4/dksi " << std::endl;

    for (int pc = 0; pc < npc; ++pc) {
        std::cout << std::setw(4) << pc + 1;
        for (int i = 0; i < 4; ++i) {
            std::cout << std::setw(15) << std::fixed << std::setprecision(6) << elem.dN_dKsi[pc][i];
        }
        std::cout << std::endl;
    }

    std::cout << "\nPochodne funkcji ksztaltu po eta:" << std::endl;
    std::cout << " pc     dN1/deta   dN2/deta   dN3/deta   dN4/deta " << std::endl;

    for (int pc = 0; pc < npc; ++pc) {
        std::cout << std::setw(4) << pc + 1;
        for (int i = 0; i < 4; ++i) {
            std::cout << std::setw(15) << std::fixed << std::setprecision(6) << elem.dN_dEta[pc][i];
        }
        std::cout << std::endl;
    }
}


void printJacobian(const Jakobian& jacobian, int npc) {

    for (int pc = 0; pc < npc; ++pc) {
        std::cout << "J_1_1: ";
        std::cout << jacobian.J[0][0][pc] << " ";
        std::cout << "J_1_2: ";
        std::cout << jacobian.J[0][1][pc] << " ";
        std::cout << "J_2_1: ";
        std::cout << jacobian.J[1][0][pc] << " ";
        std::cout << "J_2_2: ";
        std::cout << jacobian.J[1][1][pc] << " ";
        std::cout << std::endl;
    }
    for (int pc = 0; pc < npc; ++pc) {

        std::cout << "det(J) for pc " << pc + 1 << ": " << jacobian.detJ[pc] << std::endl;
        std::cout << "J_inv for pc " << pc + 1 << ":" << std::endl;
        std::cout << "J_inv_1_1: " << jacobian.J_inv[0][0][pc] << " ";
        std::cout << "J_inv_1_2: " << jacobian.J_inv[0][1][pc] << std::endl;
        std::cout << "J_inv_2_1: " << jacobian.J_inv[1][0][pc] << " ";
        std::cout << "J_inv_2_2: " << jacobian.J_inv[1][1][pc] << std::endl;
    }
}


void printCGlobalna(const std::vector<std::vector<double>>& C_globalna) {
    std::cout << "Macierz globalna: " << std::endl;
    for (size_t i = 0; i < C_globalna.size(); ++i) {  // Iteracja po wierszach
        for (size_t j = 0; j < C_globalna[i].size(); ++j) {  // Iteracja po kolumnach
            // Wypisanie elementu z okreœlon¹ precyzj¹
            std::cout << C_globalna[i][j] << " ";
        }
        std::cout << std::endl;  // Przejœcie do nowej linii po wypisaniu wiersza
    }
}


void printHbc(Element& element) {
    std::cout << "Hbc matrix for the element:" << std::endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            std::cout << element.HBC[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Vector P for the element:" << std::endl;
    for (int i = 0; i < 4; i++) {
        std::cout << element.P[i] << " ";
    }
    std::cout << std::endl;
}


void printMatrix(const  std::vector< std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double value : row) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
}


void printVector(const  std::vector<double>& vec) {
    for (double value : vec) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}
