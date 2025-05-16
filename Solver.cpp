#include <iostream>
#include <vector>
#include "Solver.h"


std::vector<double> gaussianElimination(std::vector<std::vector<double>>& H_globalna, std::vector<double>& P_global) {
    int n = H_globalna.size();
    std::vector<double> t(n); // wektor wynikoy

    // ---------------------------------------Nie dajemy na -P bo siê skompensuje wczesniej, przed calka dla [P] jest minus
    //// Przekszta³cenie wektora P na -P
    //for (int i = 0; i < n; ++i) {
    //    P_global[i] = -P_global[i];
    //}

    // Eliminacja Gaussa
    for (int k = 0; k < n; ++k) {
        // Szukanie maksymalnego elementu w kolumnie k 
        int maxIndex = k; // szukamy maksymalny element w kolumnie k
        for (int i = k + 1; i < n; ++i) {
            if (abs(H_globalna[i][k]) > abs(H_globalna[maxIndex][k])) {
                maxIndex = i;
            }
        }
        std::swap(H_globalna[k], H_globalna[maxIndex]);
        std::swap(P_global[k], P_global[maxIndex]);

        //eliminania niewiadomych ponizej glownej przekatnej 
        for (int i = k + 1; i < n; ++i) {
            double factor = H_globalna[i][k] / H_globalna[k][k];
            for (int j = k; j < n; ++j) {
                H_globalna[i][j] -= factor * H_globalna[k][j];
            }
            P_global[i] -= factor * P_global[k];
        }
    }

    // Rozwi¹zanie uk³adu za pomoc¹ podstawienia wstecznego
    for (int i = n - 1; i >= 0; --i) {
        double sum = P_global[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= H_globalna[i][j] * t[j];
        }
        t[i] = sum / H_globalna[i][i];
    }

    return t;
}


void calculateMatrix(const std::vector<std::vector<double>>& H_globalna, const std::vector<std::vector<double>>& C_globalna, double dT, std::vector<std::vector<double>>& resultMatrix) {
    // pobranie rozmiaru macierzy
    int rows = H_globalna.size();
    int cols = H_globalna[0].size();

    // Wykonanie obliczeñ: resultMatrix = H_globalna + C_globalna / dT
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            resultMatrix[i][j] = H_globalna[i][j] + C_globalna[i][j] / dT;
        }
    }
}



void calculateVector(const std::vector<double>& P_global, const std::vector<std::vector<double>>& C_globalna, const std::vector<double>& T0, double dT, std::vector<double>& resultVector) {
    int size = P_global.size();

    // Obliczenie wektora temperatury
    for (int i = 0; i < size; ++i) {
        double sum = P_global[i];

        // Dodajemy wp³yw macierzy C_globalna i temperatury T0
        for (int j = 0; j < size; ++j) {
            sum += C_globalna[i][j] * T0[j] / dT;
            // cout << "WEKTOR t0\t " << T0[j] << endl;
        }

        // Zapisujemy wynik do resultVector
        resultVector[i] = sum;
    }
}