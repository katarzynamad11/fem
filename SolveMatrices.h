#pragma once

#include <vector>
#include "Grid.h"  

// Oblicza macierz Jacobiego dla ka�dego punktu ca�kowania
void calculateJacobian(const std::vector<std::vector<double>>& dN_dKsi, const std::vector<std::vector<double>>& dN_dEta,
    const double x[4], const double y[4], Jakobian& jacobian, int npc);

// Oblicza pochodne funkcji kszta�tu po x i y
void calculateDN_dx_dy(const std::vector<std::vector<double>>& dN_dKsi, const std::vector<std::vector<double>>& dN_dEta,
    Jakobian& jacobian, int npc);

// Oblicza macierze Hpcx, Hpcy, Hpcxy dla ka�dego punktu ca�kowania
void calculate_H_matrices(const Jakobian& jacobian, H& h, int npc);

// Oblicza macierz Hpc (pomno�enie przez k oraz detJ)
void calculate_Hpc_matrix(H& h, double k, const Jakobian& jacobian, int npc);

// Oblicza lokaln� macierz H oraz sk�ada j� do globalnej H_globalna
void calculateH(H& h, const std::vector<double>& weights,int elementIndex, int N,
    Grid& grid, std::vector<std::vector<double>>& H_globalna);

// Sprawdza, czy dwa w�z�y maj� warunek brzegowy (BC)
bool hasBC(int node1, int node2, const std::vector<Node>& nodes);

// Znajduje elementy z warunkiem brzegowym (BC)
void findElementsWithBC(Grid& grid);


// Oblicza macierz C lokaln� i sk�ada do globalnej C_globalna
void calculateC(const int npc, int N, std::vector<std::vector<double>>& C_globalna,
    const Jakobian& jacobian, double c, double rho, const elementC& elC, Element& element);

// Oblicza d�ugo�� kraw�dzi pomi�dzy dwoma w�z�ami
double calculateEdgeLength(Node& n1, Node& n2);

// Oblicza macierz Hbc oraz wektor P dla elementu, a nast�pnie sk�ada do macierzy globalnych
void calculateHbc(Element& element, std::vector<Node>& nodes, double alfa, int gaussOrder, 
    double t_otoczenia, std::vector<std::vector<double>>& H_globalna, std::vector<double>& P_global);

