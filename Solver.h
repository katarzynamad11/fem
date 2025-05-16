#pragma once

//funkcja rozwi�zuj�ca uk�ad r�wna�, H globalna oraz P globalny || wynik to temperatura na ka�dym w�le 
std::vector<double> gaussianElimination(std::vector<std::vector<double>>& H_globalna, std::vector<double>& P_global);

// Oblicza macierz wynikow� [H] + [C]/dT.
void calculateMatrix(const std::vector<std::vector<double>>& H_globalna, const std::vector<std::vector<double>>& C_globalna, double dT, std::vector<std::vector<double>>& resultMatrix);

// Oblicza wektor wynikowy dla przypadku niestacjonarnego uwzgl�dniaj�c C i temperatur� pocz�tkow�.
void calculateVector(const std::vector<double>& P_global, const std::vector<std::vector<double>>& C_globalna, const std::vector<double>& T0, double dT, std::vector<double>& resultVector);