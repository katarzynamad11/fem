#include "GaussQuadrature.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// Gauss quadrature points and weights
void GaussQuadrature(int N, std::vector<double>& points, std::vector<double>& weights) {
    points.clear();
    weights.clear();

    if (N == 2) {
        points = { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };
        weights = { 1.0, 1.0 };
    }
    else if (N == 3) {
        points = { -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0) };
        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    }
    else if (N == 4) {
        points = { -0.861136, -0.339981, 0.339981, 0.861136 };
        weights = { (18.0 - std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0,
                    (18.0 + std::sqrt(30.0)) / 36.0, (18.0 - std::sqrt(30.0)) / 36.0 };
    }
    else {
        throw std::invalid_argument("Blad dla N=" + std::to_string(N));
    }
}

void calculateDerivatives(ElemUniv& elem, const std::vector<double>& ksi_points, const std::vector<double>& eta_points) {
    for (int pc = 0; pc < npc; ++pc) {
        double ksi = ksi_points[pc];
        double eta = eta_points[pc];


        elem.dN_dKsi[pc][0] = -0.25 * (1 - eta);   // dN1/dksi
        elem.dN_dKsi[pc][1] = 0.25 * (1 - eta);   // dN2/dksi
        elem.dN_dKsi[pc][2] = 0.25 * (1 + eta);   // dN3/dksi
        elem.dN_dKsi[pc][3] = -0.25 * (1 + eta);   // dN4/dksi


        elem.dN_dEta[pc][0] = -0.25 * (1 - ksi);   // dN1/deta
        elem.dN_dEta[pc][1] = -0.25 * (1 + ksi);   // dN2/deta
        elem.dN_dEta[pc][2] = 0.25 * (1 + ksi);   // dN3/deta
        elem.dN_dEta[pc][3] = 0.25 * (1 - ksi);   // dN4/deta
    }
}


void calculateShapeFunctions(std::vector<std::vector<double>>& N, int edge, int gaussOrder) {
    std::vector<double> gaussPoints, gaussWeights;
    //gaussOrder - liczba punktow calkowania = N to macierz o wymiarach [gaussOrder] x [4] bo ma 4 wêz³y
    GaussQuadrature(gaussOrder, gaussPoints, gaussWeights);

    N.resize(gaussOrder, std::vector<double>(4));

    for (int i = 0; i < gaussOrder; i++) {
        double xi = gaussPoints[i];

        switch (edge) {
        case 0: // KrawêdŸ 1: Node 1 – Node 2
            N[i][0] = 0.5 * (1 - xi); // N1
            N[i][1] = 0.5 * (1 + xi); // N2
            N[i][2] = 0.0;            // N3
            N[i][3] = 0.0;            // N4
            break;
        case 1: // KrawêdŸ 2: Node 2 – Node 3
            N[i][0] = 0.0;            // N1
            N[i][1] = 0.5 * (1 - xi); // N2
            N[i][2] = 0.5 * (1 + xi); // N3
            N[i][3] = 0.0;            // N4
            break;
        case 2: // KrawêdŸ 3: Node 3 – Node 4
            N[i][0] = 0.0;            // N1
            N[i][1] = 0.0;            // N2
            N[i][2] = 0.5 * (1 - xi); // N3
            N[i][3] = 0.5 * (1 + xi); // N4
            break;
        case 3: // KrawêdŸ 4: Node 4 – Node 1
            N[i][0] = 0.5 * (1 + xi); // N1
            N[i][1] = 0.0;            // N2
            N[i][2] = 0.0;            // N3
            N[i][3] = 0.5 * (1 - xi); // N4
            break;
        }
    }
}

void CalculateShapeFunctionMatrix(int N, elementC& elC) {
    std::vector<double> points;
    std::vector<double> weights;

    // Wywo³anie funkcji GaussQuadrature do uzyskania punktów i wag
    GaussQuadrature(N, points, weights);

    int npc = points.size(); // liczba punktów w jednym kierunku (N)

    // Liczba punktów ca³kowania (N x N)
    int totalPoints = npc * npc;

    // Alokacja pamiêci dla funkcji kszta³tu
    elC.N.resize(totalPoints, std::vector<double>(4, 0.0)); // Macierz bêdzie mia³a totalPoints wierszy, 4 kolumny

    int index = 0;
    for (double ksi : points) {
        for (double eta : points) {
            // Zapisanie wartoœci funkcji kszta³tu do odpowiednich kolumn
            elC.N[index][0] = 0.25 * (1 - ksi) * (1 - eta);  // N1
            elC.N[index][1] = 0.25 * (1 - ksi) * (1 + eta);  // N2
            elC.N[index][2] = 0.25 * (1 + ksi) * (1 + eta);  // N3
            elC.N[index][3] = 0.25 * (1 + ksi) * (1 - eta);  // N4
            index++;
        }
    }

    // Wyœwietlanie macierzy N
   /* cout << "Macierz N[" << totalPoints << "][4]:" << endl;
    for (int i = 0; i < totalPoints; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << fixed << setprecision(4) << elC.N[i][j] << " ";
        }
        cout << endl;
    }*/
}