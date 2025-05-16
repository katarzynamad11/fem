#include "SolveMatrices.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "GaussQuadrature.h"


void calculateJacobian(const std::vector<std::vector<double>>& dN_dKsi, const std::vector<std::vector<double>>& dN_dEta, const double x[4], const double y[4], Jakobian& jacobian, int npc)
{
    for (int pc = 0; pc < npc; ++pc) {
        // inicjalizacja jakobianu 2x2 wartoscia 0.0
        jacobian.J[0][0][pc] = 0.0;
        jacobian.J[0][1][pc] = 0.0;
        jacobian.J[1][0][pc] = 0.0;
        jacobian.J[1][1][pc] = 0.0;

        // jakobiany dla kazdego pc 
        for (int i = 0; i < 4; ++i) {
            jacobian.J[0][0][pc] += dN_dKsi[pc][i] * x[i];  // dN/dksi * x
            jacobian.J[0][1][pc] += dN_dKsi[pc][i] * y[i];  // dN/dksi * y
            jacobian.J[1][0][pc] += dN_dEta[pc][i] * x[i];  // dN/deta * x
            jacobian.J[1][1][pc] += dN_dEta[pc][i] * y[i];  // dN/deta * y
        }

        // dV = det J, wyznacznik
        jacobian.detJ[pc] = jacobian.J[0][0][pc] * jacobian.J[1][1][pc] - jacobian.J[0][1][pc] * jacobian.J[1][0][pc];

        // Odwrocony jakobian: 1/detJ * jacobian 
        double detJ_inv = 1.0 / jacobian.detJ[pc];
        jacobian.J_inv[0][0][pc] = jacobian.J[1][1][pc] * detJ_inv;
        jacobian.J_inv[0][1][pc] = -jacobian.J[0][1][pc] * detJ_inv;
        jacobian.J_inv[1][0][pc] = -jacobian.J[1][0][pc] * detJ_inv;
        jacobian.J_inv[1][1][pc] = jacobian.J[0][0][pc] * detJ_inv;
    }
}


void calculateDN_dx_dy(const std::vector<std::vector<double>>& dN_dKsi, const std::vector<std::vector<double>>& dN_dEta, Jakobian& jacobian, int npc) {
    //inicjalizacja dla jacobianu odwrotnego dla kazdego punktu calkowania
    for (int pc = 0; pc < npc; ++pc) {

        double invJ_00 = jacobian.J_inv[0][0][pc]; // 80
        double invJ_01 = jacobian.J_inv[0][1][pc]; // 0 
        double invJ_10 = jacobian.J_inv[1][0][pc]; // 0 
        double invJ_11 = jacobian.J_inv[1][1][pc]; // 80
        /* cout << invJ_00 << endl;
         cout << invJ_01 << endl;
         cout << invJ_10 << endl;
         cout << invJ_11 << endl;*/


         // dN/dx i dN/dy
        for (int i = 0; i < 4; ++i) {
            jacobian.dN_dX[pc][i] = dN_dKsi[pc][i] * invJ_00 + dN_dEta[pc][i] * invJ_01;
            //  cout <<"pochodna dn/dx dla pc" << pc << "\t" << jacobian.dN_dX[pc][i] << endl;
            jacobian.dN_dY[pc][i] = dN_dKsi[pc][i] * invJ_10 + dN_dEta[pc][i] * invJ_11;
            //cout << "pochodna dn/dy dla pc" << pc << "\t" << jacobian.dN_dY[pc][i] << endl;
        }
    }
}


void calculate_H_matrices(const Jakobian& jacobian, H& h, int npc) {
    // iteruje po pc 
    for (int pc = 0; pc < npc; ++pc) {
        //  Hpcx matrix
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                h.Hpcx[i][j][pc] = jacobian.dN_dX[pc][i] * jacobian.dN_dX[pc][j]; // wektor razy wektor t
                // cout << jacobian.dN_dX[pc][i] << endl;

            }
        }

        //  Hpcy matrix
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                h.Hpcy[i][j][pc] = jacobian.dN_dY[pc][i] * jacobian.dN_dY[pc][j];
            }
        }
        // hpcx + hpcy = Hpcxy (jeszcze przed wymnozeniem przed detJ i k[conductivity])
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                h.Hpcxy[i][j][pc] = h.Hpcx[i][j][pc] + h.Hpcy[i][j][pc];
            }
        }


        /* cout << "Hpcx for pc " << pc << ":" << endl;
         for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 4; ++j) {
                 cout << h.Hpcx[i][j][pc] << " ";
             }
             cout << endl;
         }

         cout << "Hpcy for pc " << pc << ":" << endl;
         for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 4; ++j) {
                 cout << h.Hpcy[i][j][pc] << " ";
             }
             cout << endl;
         }*/
    }
}

void calculate_Hpc_matrix(H& h, double k, const Jakobian& jacobian, int npc) {
    for (int pc = 0; pc < npc; ++pc) {
        double detJ = jacobian.detJ[pc];
        //  cout << "DETJ\t " << detJ << endl;
          // Mno¿enie k * detJ * Hpcxy dla ka¿dego punktu Gaussa
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                h.Hpc[i][j][pc] = k * detJ * h.Hpcxy[i][j][pc];
            }
        }


        /* cout << "Macierz Hpc dla pc " << pc << ":" << endl;
         for (int i = 0; i < 4; ++i) {
             for (int j = 0; j < 4; ++j) {
                 cout << h.Hpc[i][j][pc] << " ";
             }
             cout << endl;
         }*/
    }
}

void calculateH(H& h, const std::vector<double>& weights, int elementIndex, int N, Grid& grid, std::vector<std::vector<double>>& H_globalna) {


    // tworzymy zmienn¹ do przechowywania obliczonej macierzy H
    std::vector<std::vector<double>> H_result(4, std::vector<double>(4, 0.0));

    // iteracja przez wszystkie wartoœci npc (punkty kwadratury)
    int index = 0;  // Indeks do iteracji przez punkty kwadratury
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // Pobieramy wagi dla punktu (i, j)
            double w1 = weights[i];  // Zamiast i
            double w2 = weights[j];  // Zamiast j
            double weightProduct = w1 * w2;
            /*cout << "Hpc dla punktu " << index << endl;
            cout << "WAGAA 1 dla i = " << i << ": " << w1 << endl;
            cout << "WAGAA 2 dla j = " << j << ": " << w2 << endl;*/

            // obliczanie sk³adników macierzy Hpc (po pomno¿eniu przez odpowiednie wagi)
           /* double Hpc[4][4] = {
                {h.Hpc[0][0][index] * w1 * w2, h.Hpc[0][1][index] * w2 * w1, h.Hpc[0][2][index] * w1 * w2, h.Hpc[0][3][index] * w1 * w2},
                {h.Hpc[1][0][index] * w1 * w2, h.Hpc[1][1][index] * w2 * w1, h.Hpc[1][2][index] * w1 * w2, h.Hpc[1][3][index] * w1 * w2},
                {h.Hpc[2][0][index] * w1 * w2, h.Hpc[2][1][index] * w2 * w1, h.Hpc[2][2][index] * w1 * w2, h.Hpc[2][3][index] * w1 * w2},
                {h.Hpc[3][0][index] * w1 * w2, h.Hpc[3][1][index] * w2 * w1, h.Hpc[3][2][index] * w1 * w2, h.Hpc[3][3][index] * w1 * w2}
            };*/

            // Wyœwietlanie obliczonego Hpc dla ka¿dego punktu kwadratury
            /*cout << "Hpc dla punktu" << index << endl;
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    cout << Hpc[row][col] << " ";
                }
                cout << endl;
            }*/

            // Dodajemy do g³ównej macierzy H_result
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    // H_result[row][col] += Hpc[row][col];
                    H_result[row][col] += h.Hpc[row][col][index] * weightProduct;
                }
            }

            index++;  // Zwiêkszamy indeks dla kolejnego punktu
        }
    }

    // Wyœwietlanie koñcowej macierzy H_result
    /*cout << "Koñcowa macierz H_result:" << endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << H_result[i][j] << " ";
        }
        cout << endl;
    }*/


    h.H_matrix[elementIndex] = H_result;
    //  cout << elementIndex << endl;



    // Agregacja do macierzy globalnej - H_globalna

    const int* elementNodes = grid.El[elementIndex].ID;  // ID wêz³ów elementu
    std::cout << "ID wezlow dla elementu " << elementIndex << ": ";
    for (int i = 0; i < 4; ++i) {
        std::cout << elementNodes[i] << " ";  // Wypisujemy ka¿dy identyfikator wêz³a
    }
    std::cout << std::endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int node_i = elementNodes[i] - 1;   // Indeks wêz³a i, musimy odj¹c bo zaczynamy od 0
            int node_j = elementNodes[j] - 1;   // Indeks wêz³a j


            H_globalna[node_i][node_j] += H_result[i][j];
        }
    }

    /* cout << "Globalna macierz H_global (16x16) bez BC:" << endl;
     for (int i = 0; i < 16; ++i) {
         for (int j = 0; j < 16; ++j) {
             cout << H_globalna[i][j] << " ";
         }
         cout << endl;
     }*/


}


bool hasBC(int node1, int node2, const std::vector<Node>& nodes) {
    return (nodes[node1 - 1].BC == 1 && nodes[node2 - 1].BC == 1);
}

void findElementsWithBC(Grid& grid) {


    // iteruje po elementach i sprawdza 
    for (int e = 0; e < grid.nE; ++e) {
        Element& el = grid.El[e];

        // sprawdza wszystkie wezly, aby odnalezc te ktore s¹ na krawedziach 
        for (int i = 0; i < 4; ++i) {
            int node1 = el.ID[i];
            int node2 = el.ID[(i + 1) % 4];  // nazstepny wezel

            if (hasBC(node1, node2, grid.ND)) {
                std::cout << "Element " << e + 1 << "spelnia warunek brzegowy na: "
                    << node1 << " and " << node2 << std::endl;
            }
        }
    }
}



void calculateC( const int npc, int N, std::vector<std::vector<double>>& C_globalna, const Jakobian& jacobian, double c, double rho, const elementC& elC, Element& element) {
    // Tworzymy zmienn¹ do przechowywania macierzy Cpc dla ka¿dego punktu Gaussa
    std::vector<std::vector<std::vector<double>>> Cpc_per_point(npc, std::vector<std::vector<double>>(4, std::vector<double>(4, 0.0))); // Macierz 3D

    std::vector<double> gaussPoints, gaussWeights;
    GaussQuadrature(N, gaussPoints, gaussWeights);  // Pobranie punktów i wag Gaussa

    // Inicjalizacja macierzy C_local na 0
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            element.C_local[row][col] = 0.0;
        }
    }

    int index = 0;
    // Obliczanie Cpc dla ka¿dego punktu Gaussa
    for (int pc_x = 0; pc_x < N; ++pc_x) {
        for (int pc_y = 0; pc_y < N; ++pc_y) {
            double detJ = jacobian.detJ[index];

            // Obliczanie Cpc dla punktu Gaussa
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    // wagi Gaussa w dwóch wymiarach
                    double weight_x = gaussWeights[pc_x];
                    double weight_y = gaussWeights[pc_y];

                    // Obliczanie wartoœci Cpc
                    Cpc_per_point[index][i][j] = c * rho * (elC.N[index][i] * elC.N[index][j]) * detJ * weight_x * weight_y;
                }
            }

            // Sumowanie wartoœci Cpc do C_local
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    element.C_local[i][j] += Cpc_per_point[index][i][j];
                }
            }

            index++;
        }
    }

    // Agregacja do globalnej macierzy C (rozmiar 16x16)
    const int* elementNodes = element.ID;  // ID wêz³ów elementu

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int node_i = elementNodes[i] - 1;
            int node_j = elementNodes[j] - 1;
            C_globalna[node_i][node_j] += element.C_local[i][j];
        }
    }
}
// L to wynik obliczony z pitagorasa
double calculateEdgeLength(Node& n1, Node& n2) {
    double dx = n2.x - n1.x;
    double dy = n2.y - n1.y;
    return sqrt(dx * dx + dy * dy);
}


void calculateHbc(Element& element, std::vector<Node>& nodes, double alfa, int gaussOrder, double t_otoczenia, std::vector<std::vector<double>>& H_globalna, std::vector<double>& P_global) {
    std::vector<double> gaussPoints, gaussWeights;
    GaussQuadrature(gaussOrder, gaussPoints, gaussWeights);

    // Przechodzimy przez wszystkie krawêdzie elementu
    for (int edge = 0; edge < 4; edge++) {
        int n1 = element.ID[edge % 4] - 1; // wybiera 2 we¿³y n1 i n2 które tworz¹ dan¹ krawêdŸ 
        int n2 = element.ID[(edge + 1) % 4] - 1;

        if (nodes[n1].BC == 1 && nodes[n2].BC == 1) { // 2 wez³y musz¹ mieæ ustawion¹ flage "1"
            double Hbc_local[4][4] = { 0 };
            std::vector<double> P_local(4, 0.0);

            std::vector<std::vector<double>> N;
            calculateShapeFunctions(N, edge, gaussOrder);

            double length = calculateEdgeLength(nodes[n1], nodes[n2]);

            for (int i = 0; i < gaussOrder; i++) {
                double weight = gaussWeights[i];
                std::vector<double> Ni = N[i];

                for (int j = 0; j < 4; j++) {
                    for (int k = 0; k < 4; k++) {
                        Hbc_local[j][k] += Ni[j] * Ni[k] * weight;
                    }
                }

                for (int j = 0; j < 4; j++) {
                    P_local[j] += alfa * Ni[j] * t_otoczenia * weight * (length / 2.0);
                }
            }

            double detJ = length / 2.0;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    element.HBC[i][j] += alfa * Hbc_local[i][j] * detJ;
                }
            }

            for (int i = 0; i < 4; ++i) {
                element.P[i] += P_local[i];
            }
        }
    }

    // Dodanie HBC do globalnej macierzy H_globalna
    const int* elementNodes = element.ID;  // ID wêz³ów elementu
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int node_i = elementNodes[i] - 1;  // Indeks globalny wêz³a i
            int node_j = elementNodes[j] - 1;  // Indeks globalny wêz³a j

            H_globalna[node_i][node_j] += element.HBC[i][j];
        }
    }

    // Dodanie wektora P do globalnego wektora P_global
    for (int i = 0; i < 4; ++i) {
        int node_i = elementNodes[i] - 1;  // Indeks globalny wêz³a i
        P_global[node_i] += element.P[i];
    }

    // Wypisanie zaktualizowanej globalnej macierzy H_matrix dla elementu
    //cout << "Zaktualizowana lokalna macierz H + HBC dla elementu " << elementIndex + 1 << ":" << endl;
    //for (int i = 0; i < 4; ++i) {
    //    for (int j = 0; j < 4; ++j) {
    //        std::cout << h.H_matrix[elementIndex][i][j] << " ";
    //    }
    //    std::cout << std::endl;
    //}

    //// Wypisanie zaktualizowanej globalnej macierzy H_globalna
   /* std::cout << "Zaktualizowana globalna macierz H_globalna (16x16):" << std::endl;
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            std::cout << H_globalna[i][j] << " ";
        }
        std::cout << std::endl;
    }*/

    // Wypisanie zaktualizowanego globalnego wektora P_global
   /* cout << "Zaktualizowany globalny wektor P_global:" << endl;
    for (size_t i = 0; i < P_global.size(); ++i) {
        std::cout << P_global[i] << " ";
    }
   std::cout << std::endl;*/
}


