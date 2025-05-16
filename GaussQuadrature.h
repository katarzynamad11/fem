#pragma once
#include <vector>
#include "Grid.h"

// deklaracja funkcji kwadratury Gaussa - pomocnicza
void GaussQuadrature(int N, std::vector<double>& points, std::vector<double>& weights);

// oblicza pochodne funkcji kszta³tu
void calculateDerivatives(ElemUniv& elem, const std::vector<double>& ksi_points, const std::vector<double>& eta_points);

//oblicza funkcje ksztaltu dla BC
void calculateShapeFunctions(std::vector<std::vector<double>>& N, int edge, int gaussOrder);

//oblicza 
void CalculateShapeFunctionMatrix(int N, elementC& elC);