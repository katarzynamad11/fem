#include <iostream>
#include <vector>
#include "Grid.h"
#include "DataLoader.h"
#include "GaussQuadrature.h"
#include "SolveMatrices.h"
#include "Solver.h"
#include <algorithm>




int main() {

// Inicjalizacja struktur
Node node;
Element element;
Grid grid;
GlobalData globalData;
ElemUniv elemUniv;
// Elem2d_4 elem2d_4;
 //Surface surface;
Jakobian jacobian;
// SurfaceElement surfaceElement;
elementC elC;
H h;

// Za³adowanie elementow siatki
LoadData("Test1_4_4.txt", node, grid, globalData); // Test2_4_4_MixGrid.txt , Test1_4_4.txt , test.txt, Test3_31_31_kwadrat.txt


//Dane z pliku
int k = globalData.conductivity;
int t = globalData.tot;
int alfa = globalData.alfa;
double rho = globalData.density;
double c = globalData.specificHeat;

int N = num;
int numNodes = 4 * 4; // dla pliku testowego Test3_31_31 bedzie : 31 * 31 
// FUNKCJE KSZTALTU
std::vector<double> ksi_points, eta_points, weights;
GaussQuadrature(N, ksi_points, weights);
if (N == 2) {
    eta_points = ksi_points;

    elemUniv.ksi = { ksi_points[0], ksi_points[1], ksi_points[0], ksi_points[1] };
    elemUniv.eta = { eta_points[0], eta_points[0], eta_points[1], eta_points[1] };
    elemUniv.w1 = { weights[0], weights[1] };
}
else if (N == 3) {
    eta_points = ksi_points;

    elemUniv.ksi = {
        ksi_points[0], ksi_points[0], ksi_points[0],
        ksi_points[1], ksi_points[1], ksi_points[1],
        ksi_points[2], ksi_points[2], ksi_points[2]
    };
    elemUniv.eta = {
        eta_points[0], eta_points[1], eta_points[2],
        eta_points[0], eta_points[1], eta_points[2],
        eta_points[0], eta_points[1], eta_points[2]
    };


    elemUniv.w1 = {
   weights[0] * weights[0], weights[0] * weights[1], weights[0] * weights[2],
   weights[1] * weights[0], weights[1] * weights[1], weights[1] * weights[2],
   weights[2] * weights[0], weights[2] * weights[1], weights[2] * weights[2]
    };
}
else if (N == 4) {
    eta_points = ksi_points;


    elemUniv.ksi = {
        ksi_points[0], ksi_points[0], ksi_points[0], ksi_points[0],
        ksi_points[1], ksi_points[1], ksi_points[1], ksi_points[1],
        ksi_points[2], ksi_points[2], ksi_points[2], ksi_points[2],
        ksi_points[3], ksi_points[3], ksi_points[3], ksi_points[3]
    };
    elemUniv.eta = {
        eta_points[0], eta_points[1], eta_points[2], eta_points[3],
        eta_points[0], eta_points[1], eta_points[2], eta_points[3],
        eta_points[0], eta_points[1], eta_points[2], eta_points[3],
        eta_points[0], eta_points[1], eta_points[2], eta_points[3]
    };

    elemUniv.w1 = {
        weights[0], weights[1], weights[2], weights[3],
        weights[0], weights[1], weights[2], weights[3],
        weights[0], weights[1], weights[2], weights[3],
        weights[0], weights[1], weights[2], weights[3]
    };


}

elemUniv.dN_dKsi.resize(npc, std::vector<double>(4));  // 4 funcke ksztatltu po ksi
elemUniv.dN_dEta.resize(npc, std::vector<double>(4));  // 4 funkcje ksztaltu po eta


h.H_matrix.resize(grid.nE, std::vector<std::vector<double>>(4, std::vector<double>(4, 0.0))); // macierze lokalne  H
// vector<vector<double>> HBC(4, vector<double>(4, 0.0));



calculateDerivatives(elemUniv, elemUniv.ksi, elemUniv.eta); //oblicza pochodne funkcji ksztaltu dla macierzy H
//printDerivatives(elemUniv, npc); // wypisanie

CalculateShapeFunctionMatrix(N, elC); // oblicza funkcje ksztaltu potrzebne do obliczenia macierzy C

// Macierze i wektor globalny inicjalizujemy na 0.0
std::vector<std::vector<double>> H_globalna(numNodes, std::vector<double>(numNodes, 0.0));
std::vector<std::vector<double>>C_globalna(numNodes, std::vector<double>(numNodes, 0.0));
std::vector<double> P_global(numNodes, 0.0);


// findElementsWithBC(grid); // funckja szukajaca na ktorych wezlach/krawedziach znajduje sie warunek brzegowy

for (int el = 0; el < grid.nE; ++el) {
    std::cout << "----------------------------------------------------Dla elementu " << el + 1 << "------------------------------------------------------" << std::endl;
    Element element = grid.El[el];
    int nodeIDs[4] = { element.ID[0], element.ID[1], element.ID[2], element.ID[3] };
    double x[4], y[4];
    for (int i = 0; i < 4; ++i) {
        int nodeIndex = nodeIDs[i] - 1;
        x[i] = grid.ND[nodeIndex].x;
        y[i] = grid.ND[nodeIndex].y;
    }
    calculateJacobian(elemUniv.dN_dKsi, elemUniv.dN_dEta, x, y, jacobian, npc);
    // printJacobian(jacobian, npc);
    calculateDN_dx_dy(elemUniv.dN_dKsi, elemUniv.dN_dEta, jacobian, npc);
    calculate_H_matrices(jacobian, h, npc);
    calculate_Hpc_matrix(h, k, jacobian, npc);
    calculateH(h, weights, el, N, grid, H_globalna);
    calculateHbc(grid.El[el], grid.ND, alfa, N, t, H_globalna, P_global);
    //  printHbc(grid.El[el]);
    calculateC( npc, N, C_globalna, jacobian, c, rho, elC, element);

}



//  printCGlobalna(H_globalna);
 // printCGlobalna(C_globalna);
  // ---------------------------------------------------------------------------Rozwiazanie stacjonarne - co sie wydarzy po nieskonczenie duzym czasie-----------------------------------------------------------------------
  // warunek na zewnatrz byl 1200, dlatego caly ten obieg bedzie mial tyle 
  //cout << "********************************************** Rozwiazanie niestacjonarne **********************************************" << endl;  
  //vector<double> temperatures = gaussianElimination(H_globalna, P_global);
  //cout << "Temperatury w wezlach:" << endl;
  //for (size_t i = 0; i < temperatures.size(); ++i) {
  //    cout << "t[" << i << "] = " << temperatures[i] << endl;
  //}




  // ---------------------------------------------------------------------------Rozwiazanie niestacjonarne - obliczanie temp w odpowiednich chwilach czasowych-----------------------------------------------------------------------
  // Inicjalizacja temperatury pocz¹tkowej T0
std::vector<double> T0(numNodes, globalData.initialTemp);
double simulationStepTime = globalData.simulationStepTime;

// Macierz wynikowa [H] + [C]/dT
std::vector<std::vector<double>> resultMatrix(numNodes, std::vector<double>(numNodes, 0.0));
std::vector<double> resultVector(numNodes, 0.0);
std::vector<double> T1(numNodes, 0.0);
int numTimeSteps = static_cast<int>(globalData.simulationTime / globalData.simulationStepTime);
std::cout << "********************************************** Rozwiazanie niestacjonarne **********************************************" << std::endl;
for (int step = 0; step < numTimeSteps; ++step) {
    // obliczanie macierzy [H] + [C]/dT
    calculateMatrix(H_globalna, C_globalna, simulationStepTime, resultMatrix);

    // Wypisanie macierzy
  /*  cout << "Macierz [H] + [C]/dT dla kroku czasowego " << step + 1 << ":" << endl;
    printMatrix(resultMatrix);*/

    // obliczenie wektora {t1}
    calculateVector(P_global, C_globalna, T0, simulationStepTime, resultVector);

    //Wypisanie wektora temperatury {P}
   /* cout << "Wektor temperatury po kroku czasowym " << step + 1 << ":" << endl;
    printVector(resultVector);*/
    //Eliminacja gaussa aby obliczyc wektor x, dla konkretnego kroku czasowego
   // gaussElimination(resultMatrix, resultVector, x);
    std::vector<double> temperatures = gaussianElimination(resultMatrix, resultVector);
    // Wypisanie wyników
    /*std::cout << "Rozwi¹zanie uk³adu równañ: " << std::endl;
    for (int i = 0; i < 16; ++i) {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }*/



    // aktualizacja wektora T1
    for (int i = 0; i < numNodes; ++i) {
        // T1[i] =x[i];
        T1[i] = temperatures[i];
    }

    // zaktualizowanie T0 do nowej temperatury T1 na kolejny krok
    T0 = T1;

    // wypisanie min i max temperatury po ka¿dym kroku
    double minTemp = *min_element(temperatures.begin(), temperatures.end());
    double maxTemp = *max_element(temperatures.begin(), temperatures.end());

    /*cout << "Minimalna temperatura: " << minTemp << endl;
    cout << "Maksymalna temperatura: " << maxTemp << endl;*/

    std::cout << step << "\t" << minTemp << "\t" << maxTemp << std::endl;
}

return 0;
}