#pragma once
#include <string>
#include "Grid.h" 

// Funkcja do wypisania szczegó³ów elementów
void printElements(const Grid& grid);

void printDerivatives(const ElemUniv& elem, int npc);

void printJacobian(const Jakobian& jacobian, int npc);

void printCGlobalna(const std::vector<std::vector<double>>& C_globalna);

void printHbc(Element& element);

void printMatrix(const  std::vector< std::vector<double>>& matrix);

void printVector(const  std::vector<double>& vec);