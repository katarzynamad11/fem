#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;

void GaussQuadrature(int N, vector<double>& points, vector<double>& weights);
//------------------------------------------------------------------------------------------------STRUKTURY--------------------------------------------------------------------------
const int npc = 4;
const int num = 2;

struct Node {
    double x, y;
    int BC;  // Boundary condition
};

struct Element {
    int* ID = new int[4];
    double HBC[4][4];
    double P[4]; // wektor 
    double C_local[4][4];
};

struct Grid {
    int nN;     // liczba wezlow - 16 
    int nE;     // liczba elementow - 9 
    vector<Node> ND; // wektor przechowujący wezly
    vector<Element> El; // wektor przechowujący elementy 
};

struct GlobalData {
    double simulationTime;
    double simulationStepTime;
    double conductivity;
    double alfa;
    double tot;
    double initialTemp;
    double density;
    double specificHeat;
};

struct ElemUniv {
    vector<double> ksi;      // Wartości ksi dla punktów całkowania
    vector<double> eta;      // Wartości eta dla punktów całkowania
    vector<vector<double>> dN_dKsi; // Pochodne dN/dksi
    vector<vector<double>> dN_dEta; // Pochodne dN/deta
    vector<double> w1;      // Wartości ksi dla punktów całkowania
  

};

struct Jakobian {
    double J[2][2][npc];   
    double detJ[npc];      
    double J_inv[2][2][npc]; 

    vector<vector<double>> dN_dX; 
    vector<vector<double>> dN_dY; 

    Jakobian() {
        dN_dX.resize(npc, vector<double>(4, 0.0));  
        dN_dY.resize(npc, vector<double>(4, 0.0));  
    }
};

struct H {
    double Hpcx[4][4][npc];  // Macierz Hpcx dla każdego punktu Gaussa (4x4)
    double Hpcy[4][4][npc];  // Macierz Hpcy dla każdego punktu Gaussa (4x4)
    double Hpcxy[4][4][npc]; // Macierz Hpcxy (suma Hpcx i Hpcy)
    double Hpc[4][4][npc];   // Macierz Hpc (wynik mnożenia k * detJ * Hpcxy)
   
    double H_lok[4][4];
    double Hbc[4][4];
    vector<vector<double>> H_globalna;
    vector<vector<vector<double>>> H_matrix;
    H(int numElements, int globalSize) {
        H_matrix.resize(numElements, vector<vector<double>>(4, vector<double>(4, 0.0)));
        H_globalna.resize(globalSize, vector<double>(globalSize, 0.0));
    }
    H() {
        // Inicjalizacja macierzy 4x4 dla każdego punktu Gaussa
        for (int pc = 0; pc < npc; ++pc) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    Hpcx[i][j][pc] = 0.0;
                    Hpcy[i][j][pc] = 0.0;
                    Hpcxy[i][j][pc] = 0.0;
                    Hpc[i][j][pc] = 0.0;
                }
            }
        }
    }
};

struct elementC {
    vector<vector<double>> N;  // Macierz funkcji kształtu 
    vector<vector<double>> Cpc;  // Macierz Cpc
};


//------------------------------------------------------------------------------------------------FUNKCJE WCZYTYWANIA Z SIATKI--------------------------------------------------------------------------

void LoadData(const string& filename, Node& node, Element& element, Grid& grid, GlobalData& globalData) {
    ifstream file(filename);

    // Sprawdzanie, czy plik otworzył się poprawnie
    if (!file) {
        cout << "Błąd pliku" << endl;
        return;
    }

    double number = 0.0;
    string word;
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

    // Wczytanie danych o liczbie węzłów i elementów
    for (current_line; current_line < 10; current_line++) {
        file >> word >> word >> number;  // Pomijamy dwa słowa, a następnie wczytujemy liczbę
        data[current_line] = number;
    }

    grid.nN = static_cast<int>(data[8]); // Liczba węzłów
    grid.nE = static_cast<int>(data[9]); // Liczba elementów

    // Zmiana z alokacji dynamicznej na vector
    grid.ND.resize(grid.nN);  // Zmieniamy alokację na wektor o rozmiarze nN
    grid.El.resize(grid.nE);  // Zmieniamy alokację na wektor o rozmiarze nE

    // Wczytanie danych o węzłach
    file >> word; // Pomijamy nazwę sekcji w pliku
    for (current_line = 0; current_line < grid.nN; current_line++) {
        file >> number >> word >> number >> word;  // Pomijamy dane, a następnie wczytujemy współrzędne
        grid.ND[current_line].x = number;
        file >> number;
        grid.ND[current_line].y = number;
    }

    // Wczytanie danych o elementach
    file >> word >> word; // Pomijamy nagłówek
    for (current_line = 0; current_line < grid.nE; current_line++) {
        file >> number >> word >> number >> word; // Pomijamy dane, a następnie wczytujemy ID węzłów elementu
        grid.El[current_line].ID[0] = static_cast<int>(number); // jawnie przekonwertuje na typ int 
        file >> number >> word;
        grid.El[current_line].ID[1] = static_cast<int>(number);
        file >> number >> word;
        grid.El[current_line].ID[2] = static_cast<int>(number);
        file >> number;
        grid.El[current_line].ID[3] = static_cast<int>(number);
    }

    // Wczytanie warunków brzegowych (BC) dla węzłów
    file >> word;
    int index;
    while (file >> index) {
        // Jeżeli istnieje wartość po indeksie, czytamy ją
        if (file.peek() != '\n') {  // sprawdzamy, czy po numerze jest jeszcze tekst
            file >> word;
        }
        else {
            word = "";  // jeżeli brak word, przypisujemy pustą wartość
        }

        //cout << "Wczytano index: " << index << ", word: " << word << endl;

        // Wczytujemy ID węzła i przypisujemy warunek brzegowy
        if (index > 0 && index <= grid.nN) {
            grid.ND[index - 1].BC = 1; // Indeksowanie od 1, więc odejmujemy 1
           // cout << "Przypisano BC=1 do węzła " << index << endl;
        }
        else {
            cout << "Błędny numer węzła w warunkach brzegowych: " << index << endl;
        }
    }



    // Przypisanie pierwszego węzła do przekazanego argumentu 'node'
    if (grid.nN > 1) {
        node = grid.ND[1]; // Odczytujemy drugi węzeł (indeksowanie od 1)
    }
    else {
        cout << "Brak wystarczającej liczby węzłów!" << endl;
    }

    file.close();  // Zamknięcie pliku
}

// Funkcja do wypisania szczegółów elementów
void printElements(const Grid& grid) {
    for (int i = 0; i < grid.nE; ++i) {
        cout << "Element ID: " << i + 1 << endl;
        cout << "Nodes: ";
        for (int j = 0; j < 4; ++j) {
            int node_id = grid.El[i].ID[j] - 1; // Indeksowanie węzłów od 1 w pliku, więc odejmujemy 1
            const Node& node = grid.ND[node_id];
            cout << "(" << node.x << ", " << node.y << ") ";
        }
        cout << endl;
    }
    for (int i = 0; i < grid.ND.size(); ++i) {
        const Node& node = grid.ND[i];  // Dostęp do konkretnego węzła
        cout << "Node ID: " << i + 1 << " BC: " << node.BC << endl;
    }
}

//------------------------------------------------------------------------------------------------FUNKCJE OBLICZANIA POCHODNYCH KSZTALTU--------------------------------------------------------------------------

void GaussQuadrature(int N, vector<double>& points, vector<double>& weights) {
    points.clear();
    weights.clear();

    if (N == 2) { // 2-punktowy Gauss quadrature
        points = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
        weights = { 1.0, 1.0 , 1.0, 1.0 };
    }
    else if (N == 3) { // 3-punktowy Gauss quadrature
        points = { -sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0) };
        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    }
    else if (N == 4) { // 4-punktowy Gauss quadrature
        points = { -0.861136, -0.339981, 0.339981, 0.861136 };
        weights = { (18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0,
                    (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0 };
    }
    else {
        throw invalid_argument("Blad dla N=" + to_string(N));
    }
}

void calculateDerivatives(ElemUniv& elem, const vector<double>& ksi_points, const vector<double>& eta_points) {
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

void printDerivatives(const ElemUniv& elem, int npc) {
    
 
    cout << "Integration points (N=" << npc << "):" << endl;
    for (size_t i = 0; i < elem.ksi.size(); ++i) {
        cout << "Point " << i + 1 << ": ksi = " << elem.ksi[i] << ", eta = " << elem.eta[i] << endl;
    }

    cout << "\nPochodne funkcji ksztaltu po ksi:" << endl;
    cout << " pc     dN1/dksi    dN2/dksi    dN3/dksi    dN4/dksi " << endl;

    for (int pc = 0; pc < npc; ++pc) {
        cout << setw(4) << pc + 1;
        for (int i = 0; i < 4; ++i) {
            cout << setw(15) << fixed << setprecision(6) << elem.dN_dKsi[pc][i];
        }
        cout << endl;
    }

    cout << "\nPochodne funkcji ksztaltu po eta:" << endl;
    cout << " pc     dN1/deta   dN2/deta   dN3/deta   dN4/deta " << endl;

    for (int pc = 0; pc < npc; ++pc) {
        cout << setw(4) << pc + 1;
        for (int i = 0; i < 4; ++i) {
            cout << setw(15) << fixed << setprecision(6) << elem.dN_dEta[pc][i];
        }
        cout << endl;
    }
}

//---------------------------------------------------------------------------------------------------FUNKCJE OBLICZAJACE JAKOBIANY----------------------------------------------------------------------------------

//oblicza macierz Jacobiego dla każdego punktu całkowania
void calculateJacobian(const vector<vector<double>>& dN_dKsi, const vector<vector<double>>& dN_dEta,  const double x[4], const double y[4], Jakobian& jacobian, int npc)
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

void printJacobian(const Jakobian& jacobian, int npc) {

    for (int pc = 0; pc < npc; ++pc) {
        cout << "J_1_1: ";
        cout << jacobian.J[0][0][pc] << " ";
        cout << "J_1_2: ";
        cout << jacobian.J[0][1][pc] << " ";
        cout << "J_2_1: ";
        cout << jacobian.J[1][0][pc] << " ";
        cout << "J_2_2: ";
        cout << jacobian.J[1][1][pc] << " ";
        cout << endl;
    }
    for (int pc = 0; pc < npc; ++pc) {

        cout << "det(J) for pc " << pc + 1 << ": " << jacobian.detJ[pc] << endl;
        cout << "J_inv for pc " << pc + 1 << ":" << endl;
        cout << "J_inv_1_1: " << jacobian.J_inv[0][0][pc] << " ";
        cout << "J_inv_1_2: " << jacobian.J_inv[0][1][pc] << endl;
        cout << "J_inv_2_1: " << jacobian.J_inv[1][0][pc] << " ";
        cout << "J_inv_2_2: " << jacobian.J_inv[1][1][pc] << endl;
    }
}

//---------------------------------------------------------------------------------------------------FUNKCJE OBLICZAJACE DN/DX i DN/DY----------------------------------------------------------------------------------

void calculateDN_dx_dy(const vector<vector<double>>& dN_dKsi, const vector<vector<double>>& dN_dEta, Jakobian& jacobian, int npc) {
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

//---------------------------------------------------------------------------------------------------HPC i H globalna----------------------------------------------------------------------------------


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
          // Mnożenie k * detJ * Hpcxy dla każdego punktu Gaussa
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

void calculateH(H& h, const vector<double>& weights, const int npc, int elementIndex, int N, Grid& grid, vector<vector<double>>& H_globalna) {


    // tworzymy zmienną do przechowywania obliczonej macierzy H
    vector<vector<double>> H_result(4, vector<double>(4, 0.0));

    // iteracja przez wszystkie wartości npc (punkty kwadratury)
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

            // obliczanie składników macierzy Hpc (po pomnożeniu przez odpowiednie wagi)
           /* double Hpc[4][4] = {
                {h.Hpc[0][0][index] * w1 * w2, h.Hpc[0][1][index] * w2 * w1, h.Hpc[0][2][index] * w1 * w2, h.Hpc[0][3][index] * w1 * w2},
                {h.Hpc[1][0][index] * w1 * w2, h.Hpc[1][1][index] * w2 * w1, h.Hpc[1][2][index] * w1 * w2, h.Hpc[1][3][index] * w1 * w2},
                {h.Hpc[2][0][index] * w1 * w2, h.Hpc[2][1][index] * w2 * w1, h.Hpc[2][2][index] * w1 * w2, h.Hpc[2][3][index] * w1 * w2},
                {h.Hpc[3][0][index] * w1 * w2, h.Hpc[3][1][index] * w2 * w1, h.Hpc[3][2][index] * w1 * w2, h.Hpc[3][3][index] * w1 * w2}
            };*/

            // Wyświetlanie obliczonego Hpc dla każdego punktu kwadratury
            /*cout << "Hpc dla punktu" << index << endl;
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                    cout << Hpc[row][col] << " ";
                }
                cout << endl;
            }*/

            // Dodajemy do głównej macierzy H_result
            for (int row = 0; row < 4; ++row) {
                for (int col = 0; col < 4; ++col) {
                   // H_result[row][col] += Hpc[row][col];
                    H_result[row][col] += h.Hpc[row][col][index] * weightProduct;
                }
            }

            index++;  // Zwiększamy indeks dla kolejnego punktu
        }
    }

    // Wyświetlanie końcowej macierzy H_result
    /*cout << "Końcowa macierz H_result:" << endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << H_result[i][j] << " ";
        }
        cout << endl;
    }*/

    
    h.H_matrix[elementIndex] = H_result;
    //  cout << elementIndex << endl;



    // Agregacja do macierzy globalnej - H_globalna

    const int* elementNodes = grid.El[elementIndex].ID;  // ID węzłów elementu
    cout << "ID węzłów dla elementu " << elementIndex << ": ";
    for (int i = 0; i < 4; ++i) {
        cout << elementNodes[i] << " ";  // Wypisujemy każdy identyfikator węzła
    }
    cout << endl;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int node_i = elementNodes[i] - 1;   // Indeks węzła i, musimy odjąc bo zaczynamy od 0
            int node_j = elementNodes[j] - 1;   // Indeks węzła j


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


//---------------------------------------------------------------------------------------------------HBC----------------------------------------------------------------------------------

bool hasBC(int node1, int node2, const vector<Node>& nodes) {
    return (nodes[node1 - 1].BC == 1 && nodes[node2 - 1].BC == 1); 
}

void findElementsWithBC(Grid& grid) {


    // iteruje po elementach i sprawdza 
    for (int e = 0; e < grid.nE; ++e) {
        Element& el = grid.El[e];

        // sprawdza wszystkie wezly, aby odnalezc te ktore są na krawedziach 
        for (int i = 0; i < 4; ++i) {
            int node1 = el.ID[i];
            int node2 = el.ID[(i + 1) % 4];  // nazstepny wezel

            if (hasBC(node1, node2, grid.ND)) {
                cout << "Element " << e + 1 <<  "spelnia warunek brzegowy na: "
                    << node1 << " and " << node2 << endl;
            }
        }
    }
}


//void calculateShapeFunctions(double ksi, double eta, double N[4]) {
//    N[0] = 0.25 * (1 - ksi) * (1 - eta);  // N1
//    N[1] = 0.25 * (1 + ksi) * (1 - eta);  // N2
//    N[2] = 0.25 * (1 + ksi) * (1 + eta);  // N3
//    N[3] = 0.25 * (1 - ksi) * (1 + eta);  // N4
//}
//----------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------FUNCKJE KSZTALTU N dla HBC-------------------------------------------------------------------------------------
//oblicza funkcje ksztaltu dla BC
void calculateShapeFunctions(vector<vector<double>>& N, int edge, int gaussOrder) {
    vector<double> gaussPoints, gaussWeights;
    //gaussOrder - liczba punktow calkowania = N to macierz o wymiarach [gaussOrder] x [4] bo ma 4 węzły
    GaussQuadrature(gaussOrder, gaussPoints, gaussWeights);

    N.resize(gaussOrder, vector<double>(4)); 

    for (int i = 0; i < gaussOrder; i++) {
        double xi = gaussPoints[i];

        switch (edge) {
        case 0: // Krawędź 1: Node 1 – Node 2
            N[i][0] = 0.5 * (1 - xi); // N1
            N[i][1] = 0.5 * (1 + xi); // N2
            N[i][2] = 0.0;            // N3
            N[i][3] = 0.0;            // N4
            break;
        case 1: // Krawędź 2: Node 2 – Node 3
            N[i][0] = 0.0;            // N1
            N[i][1] = 0.5 * (1 - xi); // N2
            N[i][2] = 0.5 * (1 + xi); // N3
            N[i][3] = 0.0;            // N4
            break;
        case 2: // Krawędź 3: Node 3 – Node 4
            N[i][0] = 0.0;            // N1
            N[i][1] = 0.0;            // N2
            N[i][2] = 0.5 * (1 - xi); // N3
            N[i][3] = 0.5 * (1 + xi); // N4
            break;
        case 3: // Krawędź 4: Node 4 – Node 1
            N[i][0] = 0.5 * (1 + xi); // N1
            N[i][1] = 0.0;            // N2
            N[i][2] = 0.0;            // N3
            N[i][3] = 0.5 * (1 - xi); // N4
            break;
        }
    }
}


//------------------------------------------------------------------------------------------FUNCKJE KSZTALTU N dla macierzy C oraz macierz C lokalna i globalna -------------------------------------------------------------------------------------

void CalculateShapeFunctionMatrix(int N, elementC& elC) {
    vector<double> points;
    vector<double> weights;

    // Wywołanie funkcji GaussQuadrature do uzyskania punktów i wag
    GaussQuadrature(N, points, weights);

    int npc = points.size(); // liczba punktów w jednym kierunku (N)

    // Liczba punktów całkowania (N x N)
    int totalPoints = npc * npc;

    // Alokacja pamięci dla funkcji kształtu
    elC.N.resize(totalPoints, vector<double>(4, 0.0)); // Macierz będzie miała totalPoints wierszy, 4 kolumny

    int index = 0;
    for (double ksi : points) {
        for (double eta : points) {
            // Zapisanie wartości funkcji kształtu do odpowiednich kolumn
            elC.N[index][0] = 0.25 * (1 - ksi) * (1 - eta);  // N1
            elC.N[index][1] = 0.25 * (1 - ksi) * (1 + eta);  // N2
            elC.N[index][2] = 0.25 * (1 + ksi) * (1 + eta);  // N3
            elC.N[index][3] = 0.25 * (1 + ksi) * (1 - eta);  // N4
            index++;
        }
    }

    // Wyświetlanie macierzy N
   /* cout << "Macierz N[" << totalPoints << "][4]:" << endl;
    for (int i = 0; i < totalPoints; ++i) {
        for (int j = 0; j < 4; ++j) {
            cout << fixed << setprecision(4) << elC.N[i][j] << " ";
        }
        cout << endl;
    }*/
}



void calculateC(H& h, const int npc, int N, Grid& grid, vector<vector<double>>& C_globalna, const Jakobian& jacobian, double c, double rho, const elementC& elC, Element& element) {
    // Tworzymy zmienną do przechowywania macierzy Cpc dla każdego punktu Gaussa
    vector<vector<vector<double>>> Cpc_per_point(npc, vector<vector<double>>(4, vector<double>(4, 0.0))); // Macierz 3D

    vector<double> gaussPoints, gaussWeights;
    GaussQuadrature(N, gaussPoints, gaussWeights);  // Pobranie punktów i wag Gaussa

    // Inicjalizacja macierzy C_local na 0
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            element.C_local[row][col] = 0.0; 
        }
    }

    int index = 0;
    // Obliczanie Cpc dla każdego punktu Gaussa
    for (int pc_x = 0; pc_x < N; ++pc_x) {
        for (int pc_y = 0; pc_y < N; ++pc_y) {
            double detJ = jacobian.detJ[index];  

            // Obliczanie Cpc dla punktu Gaussa
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    // wagi Gaussa w dwóch wymiarach
                    double weight_x = gaussWeights[pc_x];
                    double weight_y = gaussWeights[pc_y];

                    // Obliczanie wartości Cpc
                    Cpc_per_point[index][i][j] = c * rho * (elC.N[index][i] * elC.N[index][j]) * detJ * weight_x * weight_y;
                }
            }

            // Sumowanie wartości Cpc do C_local
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    element.C_local[i][j] += Cpc_per_point[index][i][j];
                }
            }

            index++; 
        }
    }

    // Agregacja do globalnej macierzy C (rozmiar 16x16)
    const int* elementNodes = element.ID;  // ID węzłów elementu

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int node_i = elementNodes[i] - 1;
            int node_j = elementNodes[j] - 1;
            C_globalna[node_i][node_j] += element.C_local[i][j];
        }
    }
}

void printCGlobalna(const std::vector<std::vector<double>>& C_globalna) {
    cout << "Macierz globalna: " << endl;
    for (int i = 0; i < C_globalna.size(); ++i) {  // Iteracja po wierszach
        for (int j = 0; j < C_globalna[i].size(); ++j) {  // Iteracja po kolumnach
            // Wypisanie elementu z określoną precyzją
            std::cout << C_globalna[i][j] << " ";
        }
        std::cout << std::endl;  // Przejście do nowej linii po wypisaniu wiersza
    }
}


//-------------------------------------------------------------------------------------------- dlugosc boku: L/2  ----------------------------------------------------------------------------------------------------
// L to wynik obliczony z pitagorasa
double calculateEdgeLength(Node& n1, Node& n2) {
    double dx = n2.x - n1.x;
    double dy = n2.y - n1.y;
    return sqrt(dx * dx + dy * dy);
}


//------------------------------------------------------------------------------------Macierz HBC, wektor P---------------------------------------------------------------------------------------------------------
void calculateHbc(Element& element, vector<Node>& nodes, double alfa, int gaussOrder, H& h, double t_otoczenia, int elementIndex, vector<vector<double>>& H_globalna, vector<double>& P_global) {
    vector<double> gaussPoints, gaussWeights;
    GaussQuadrature(gaussOrder, gaussPoints, gaussWeights); 

    // Przechodzimy przez wszystkie krawędzie elementu
    for (int edge = 0; edge < 4; edge++) {
        int n1 = element.ID[edge % 4] - 1; // wybiera 2 weżły n1 i n2 które tworzą daną krawędź 
        int n2 = element.ID[(edge + 1) % 4] - 1;

        if (nodes[n1].BC == 1 && nodes[n2].BC == 1) { // 2 wezły muszą mieć ustawioną flage "1"
            double Hbc_local[4][4] = { 0 };
            vector<double> P_local(4, 0.0);

            vector<vector<double>> N;
            calculateShapeFunctions(N, edge, gaussOrder);
           
            double length = calculateEdgeLength(nodes[n1], nodes[n2]);

            for (int i = 0; i < gaussOrder; i++) {
                double weight = gaussWeights[i];
                vector<double> Ni = N[i];

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
    const int* elementNodes = element.ID;  // ID węzłów elementu
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int node_i = elementNodes[i] - 1;  // Indeks globalny węzła i
            int node_j = elementNodes[j] - 1;  // Indeks globalny węzła j

            H_globalna[node_i][node_j] += element.HBC[i][j];
        }
    }

    // Dodanie wektora P do globalnego wektora P_global
    for (int i = 0; i < 4; ++i) {
        int node_i = elementNodes[i] - 1;  // Indeks globalny węzła i
        P_global[node_i] += element.P[i];
    }

    // Wypisanie zaktualizowanej globalnej macierzy H_matrix dla elementu
    //cout << "Zaktualizowana lokalna macierz H + HBC dla elementu " << elementIndex + 1 << ":" << endl;
    //for (int i = 0; i < 4; ++i) {
    //    for (int j = 0; j < 4; ++j) {
    //        cout << h.H_matrix[elementIndex][i][j] << " ";
    //    }
    //    cout << endl;
    //}

    //// Wypisanie zaktualizowanej globalnej macierzy H_globalna
   /* cout << "Zaktualizowana globalna macierz H_globalna (16x16):" << endl;
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            cout << H_globalna[i][j] << " ";
        }
        cout << endl;
    }*/

    // Wypisanie zaktualizowanego globalnego wektora P_global
   /* cout << "Zaktualizowany globalny wektor P_global:" << endl;
    for (size_t i = 0; i < P_global.size(); ++i) {
        cout << P_global[i] << " ";
    }
    cout << endl;*/
}


// Funkcja do wypisania macierzy Hbc elementu
void printHbc(Element& element) {
    cout << "Hbc matrix for the element:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << element.HBC[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Vector P for the element:" << endl;
    for (int i = 0; i < 4; i++) {
        cout << element.P[i] << " ";
    }
    cout << endl;
}


//------------------------------------------------------------------------------------------------- [H]{t} + [P] = 0 -------------------------------------------------------------------------------------------
//funkcja rozwiązująca układ równań, H globalna oraz P globalny || wynik to temperatura na każdym węźle 
vector<double> gaussianElimination(vector<vector<double>>& H_globalna, vector<double>& P_global) {
    int n = H_globalna.size();
    vector<double> t(n); // wektor wynikoy

    // ---------------------------------------Nie dajemy na -P bo się skompensuje wczesniej, przed calka dla [P] jest minus
    //// Przekształcenie wektora P na -P
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
        swap(H_globalna[k], H_globalna[maxIndex]);
        swap(P_global[k], P_global[maxIndex]);

        //eliminania niewiadomych ponizej glownej przekatnej 
        for (int i = k + 1; i < n; ++i) { 
            double factor = H_globalna[i][k] / H_globalna[k][k];
            for (int j = k; j < n; ++j) {
                H_globalna[i][j] -= factor * H_globalna[k][j];
            }
            P_global[i] -= factor * P_global[k];
        }
    }

    // Rozwiązanie układu za pomocą podstawienia wstecznego
    for (int i = n - 1; i >= 0; --i) {
        double sum = P_global[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= H_globalna[i][j] * t[j];
        }
        t[i] = sum / H_globalna[i][i];
    }

    return t;
}


//------------------------------------------------------------------------------------------------- Rozwiazanie niestacjonarne  -------------------------------------------------------------------------------------------

void calculateMatrix(const vector<vector<double>>& H_globalna,  const vector<vector<double>>& C_globalna,  double dT, vector<vector<double>>& resultMatrix) {
    // pobranie rozmiaru macierzy
    int rows = H_globalna.size();
    int cols = H_globalna[0].size();

    // Wykonanie obliczeń: resultMatrix = H_globalna + C_globalna / dT
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            resultMatrix[i][j] = H_globalna[i][j] + C_globalna[i][j] / dT;
        }
    }
}

void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double value : row) {
            cout << value << " ";
        }
        cout << endl;
    }
}

void calculateVector(const vector<double>& P_global, const vector<vector<double>>& C_globalna, const vector<double>& T0, double dT, vector<double>& resultVector) {
    int size = P_global.size();

    // Obliczenie wektora temperatury
    for (int i = 0; i < size; ++i) {
        double sum = P_global[i];

        // Dodajemy wpływ macierzy C_globalna i temperatury T0
        for (int j = 0; j < size; ++j) {
            sum += C_globalna[i][j] * T0[j] / dT;
           // cout << "WEKTOR t0\t " << T0[j] << endl;
        }

        // Zapisujemy wynik do resultVector
        resultVector[i] = sum;
    }
}

void printVector(const vector<double>& vec) {
    for (double value : vec) {
        cout << value << " ";
    }
    cout << endl;
}






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

    // Załadowanie elementow siatki
    LoadData("Test2_4_4_MixGrid.txt", node, element, grid, globalData); // Test2_4_4_MixGrid.txt , Test1_4_4.txt , test.txt, Test3_31_31_kwadrat.txt


    //Dane z pliku
    int k = globalData.conductivity;
    int t = globalData.tot;
    int alfa = globalData.alfa;
    double rho = globalData.density;
    double c = globalData.specificHeat;

    int N = num; 
    int numNodes = 4 * 4; // dla pliku testowego Test3_31_31 bedzie : 31 * 31 
    // FUNKCJE KSZTALTU
    vector<double> ksi_points, eta_points, weights;
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

   elemUniv.dN_dKsi.resize(npc, vector<double>(4));  // 4 funcke ksztatltu po ksi
   elemUniv.dN_dEta.resize(npc, vector<double>(4));  // 4 funkcje ksztaltu po eta
   
  
    h.H_matrix.resize(grid.nE, vector<vector<double>>(4, vector<double>(4, 0.0))); // macierze lokalne  H
   // vector<vector<double>> HBC(4, vector<double>(4, 0.0));

   
    
    calculateDerivatives(elemUniv, elemUniv.ksi, elemUniv.eta); //oblicza pochodne funkcji ksztaltu dla macierzy H
    //printDerivatives(elemUniv, npc); // wypisanie

    CalculateShapeFunctionMatrix(N, elC); // oblicza funkcje ksztaltu potrzebne do obliczenia macierzy C
   
    // Macierze i wektor globalny inicjalizujemy na 0.0
    vector<vector<double>> H_globalna(numNodes, vector<double>(numNodes, 0.0));
    vector<vector<double>>C_globalna(numNodes, vector<double>(numNodes, 0.0));
    vector<double> P_global(numNodes, 0.0);

   
    // findElementsWithBC(grid); // funckja szukajaca na ktorych wezlach/krawedziach znajduje sie warunek brzegowy

    for (int el = 0; el < grid.nE; ++el) {
        cout << "----------------------------------------------------Dla elementu " << el + 1 << "------------------------------------------------------" << endl;
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
        calculateH(h, weights, npc, el, N, grid, H_globalna);
        calculateHbc(grid.El[el], grid.ND, alfa,N, h,t ,el, H_globalna, P_global);
      //  printHbc(grid.El[el]);
        calculateC(h, npc, N, grid, C_globalna, jacobian, c, rho, elC, element);

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
    // Inicjalizacja temperatury początkowej T0
    vector<double> T0(numNodes, globalData.initialTemp);
    double simulationStepTime = globalData.simulationStepTime;

    // Macierz wynikowa [H] + [C]/dT
    vector<vector<double>> resultMatrix(numNodes, vector<double>(numNodes, 0.0));
    vector<double> resultVector(numNodes, 0.0);
    vector<double> T1(numNodes, 0.0);
    int numTimeSteps = static_cast<int>(globalData.simulationTime / globalData.simulationStepTime);
    cout << "********************************************** Rozwiazanie niestacjonarne **********************************************" << endl;
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
        vector<double> temperatures = gaussianElimination(resultMatrix, resultVector);
        // Wypisanie wyników
        /*std::cout << "Rozwiązanie układu równań: " << std::endl;
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

        // wypisanie min i max temperatury po każdym kroku
        double minTemp = *min_element(temperatures.begin(), temperatures.end());
        double maxTemp = *max_element(temperatures.begin(), temperatures.end());

        /*cout << "Minimalna temperatura: " << minTemp << endl;
        cout << "Maksymalna temperatura: " << maxTemp << endl;*/

        cout << step << "\t" << minTemp << "\t" << maxTemp << endl;
    }

    return 0;
}
