#include <iostream>
#include "DataLoader.h"
#include <fstream>
#include <string>

void LoadData(const std::string& filename, Node& node, Grid& grid, GlobalData& globalData) {
    std::ifstream file(filename);

    // Sprawdzanie, czy plik otworzy³ siê poprawnie
    if (!file) {
        std::cout << "B³¹d pliku" << std::endl;
        return;
    }

    double number = 0.0;
    std::string word;
    double data[10];  // Tablica pomocnicza do przechowywania danych

    int current_line = 0;

    // Wczytanie pierwszych 8 parametrów do danych globalnych
    for (current_line = 0; current_line < 8; current_line++) {
        file >> word >> number;
        data[current_line] = number;
    }

    globalData.simulationTime = data[0];
    globalData.simulationStepTime = data[1];
    globalData.conductivity = data[2];
    globalData.alfa = data[3];
    globalData.tot = data[4];
    globalData.initialTemp = data[5];
    globalData.density = data[6];
    globalData.specificHeat = data[7];

    // Wczytanie danych o liczbie wêz³ów i elementów
    for (; current_line < 10; current_line++){
        file >> word >> word >> number;  // Pomijamy dwa s³owa, a nastêpnie wczytujemy liczbê
        data[current_line] = number;
    }

    grid.nN = static_cast<int>(data[8]); // Liczba wêz³ów
    grid.nE = static_cast<int>(data[9]); // Liczba elementów

    // Zmiana z alokacji dynamicznej na vector
    grid.ND.resize(grid.nN);  // Zmieniamy alokacjê na wektor o rozmiarze nN
    grid.El.resize(grid.nE);  // Zmieniamy alokacjê na wektor o rozmiarze nE

    // Wczytanie danych o wêz³ach
    file >> word; // Pomijamy nazwê sekcji w pliku
    for (current_line = 0; current_line < grid.nN; current_line++) {
        file >> number >> word >> number >> word;  // Pomijamy dane, a nastêpnie wczytujemy wspó³rzêdne
        grid.ND[current_line].x = number;
        file >> number;
        grid.ND[current_line].y = number;
    }

    // Wczytanie danych o elementach
    file >> word >> word; // Pomijamy nag³ówek
    for (current_line = 0; current_line < grid.nE; current_line++) {
        file >> number >> word >> number >> word; // Pomijamy dane, a nastêpnie wczytujemy ID wêz³ów elementu
        grid.El[current_line].ID[0] = static_cast<int>(number); // jawnie przekonwertuje na typ int 
        file >> number >> word;
        grid.El[current_line].ID[1] = static_cast<int>(number);
        file >> number >> word;
        grid.El[current_line].ID[2] = static_cast<int>(number);
        file >> number;
        grid.El[current_line].ID[3] = static_cast<int>(number);
    }

    // Wczytanie warunków brzegowych (BC) dla wêz³ów
    file >> word;
    int index;
    while (file >> index) {
        // Je¿eli istnieje wartoœæ po indeksie, czytamy j¹
        if (file.peek() != '\n') {  // sprawdzamy, czy po numerze jest jeszcze tekst
            file >> word;
        }
        else {
            word = "";  // je¿eli brak word, przypisujemy pust¹ wartoœæ
        }

        //cout << "Wczytano index: " << index << ", word: " << word << endl;

        // Wczytujemy ID wêz³a i przypisujemy warunek brzegowy
        if (index > 0 && index <= grid.nN) {
            grid.ND[index - 1].BC = 1; // Indeksowanie od 1, wiêc odejmujemy 1
            // cout << "Przypisano BC=1 do wêz³a " << index << endl;
        }
        else {
            std::cout << "B³êdny numer wêz³a w warunkach brzegowych: " << index << std::endl;
        }
    }



    // Przypisanie pierwszego wêz³a do przekazanego argumentu 'node'
    if (grid.nN > 1) {
        node = grid.ND[1]; // Odczytujemy drugi wêze³ (indeksowanie od 1)
    }
    else {
        std::cout << "Brak wystarczaj¹cej liczby wêz³ów!" << std::endl;
    }

    file.close();  // Zamkniêcie pliku
}