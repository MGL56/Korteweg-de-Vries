#include "maclasse.h"
#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>
#include <vector>
#include <chrono>

int main() {

    // Paramètres du problème
    int N = 100; // nombre d'intervalles
    double dt = 0.01; // pas de temps
    int n_time = 1500; // Nombre de pas de temps
    double eps = 0.1; 
    double L = 10.0; // Demi-longueur de l'espace
    double dx = 2*L/N; // pas d'espace


    // Initialisation des vecteurs X_n+1, X_n, X_n-1 et du maillage
    std::vector <double> xi_2, xi_1, xi_0, x_space;

    for (int i = 0; i < N; i++){
        x_space.push_back(-L + 2*(double)i*L/N); // Création du maillage xi = -L + 2*i*L/N
        xi_2.push_back(0); // Initialisation de X_n+1 = 0
        xi_1.push_back(exp(-2*pow(x_space[i], 2))); // Condition initialle : exp(-2x²)
        xi_0.push_back(exp(-2*pow(x_space[i], 2))); // X_-1 = X_0
    }

    // On cherche à résoudre l'équation (I + dt/2M)X_n+1 = (I-dt/2*M)X_n
    // avec M = D1 + eps/6*D3 + eps/2*(diag(z_n+1/2)*D1 + D1*diag(z_n+1/2))
    // où diag(z_n+1/2) est une matrice diagonale dont la diagonale est égale
    // à z_n+1/2
    // Avec z_n+1/2 = 1/2*(3X_n - X_n-1)


    // Création de la matrice Identité    
    MaClasse I; 
    I.Reallocate(N); // Allocation des tableaux
    I.Zero(); //  Initialise les valeurs nont nulles à 0
    for (int i=0; i<N; i=i+1)
        I.Set(i, i, 1); // Rempli la diagonale de la valeur 1

    // Création de la matrice D1
    MaClasse D1;
    D1.Reallocate(N);
    D1.Zero();
    D1.Set(0, N-1, -1);
    D1.Set(N-1, 0, 1);
    for (int i=0; i<N; i++){
        D1.Set(i-1, i, 1);
        D1.Set(i, i-1, -1);

    }
    double coeff_D1 = (1.0/2)*(1.0/dx);
    MaClasse D1_ = D1*coeff_D1;

    // Création de la matrice D3
    MaClasse D3;
    D3.Reallocate(N);
    D3.Zero();
    D3.Set(0, N-1, 2);
    D3.Set(N-1, 0, -2);
    D3.Set(0, N-2, -1);
    D3.Set(1, N-1, -1);
    D3.Set(N-2, 0, 1);
    D3.Set(N-1, 1, 1);
    
    for (int i=0; i<N-1; i=i+1){
        D3.Set(i, i+1, -2);
        D3.Set(i+1, i, 2);
        if (i<N-2){
            D3.Set(i, i+2, 1);
            D3.Set(i+2, i, -1);
        }
    }
    double coeff_D3 = (1.0/2)*(1.0/dx)*(1.0/dx)*(1.0/dx);
    MaClasse D3_ = D3*coeff_D3;


    // Matrice de stockage des résultats
    std::vector<std::vector<double>> results(n_time + 1 , vector<double> (N, 0));
    for (int k=0; k<x_space.size(); k++)
        results[0][k] = xi_1[k];
    

    // Résolution du problème par itération succéssive
    for (int i=0; i<n_time; i++){
        std::vector <double> z(N, 0);
        for(int i=0; i<N; i=i+1)
            z[i] = 1/2*(3*xi_1[i] - xi_0[i]); // Création du vecteur Z_n
        
        MaClasse facteur;
        facteur.Reallocate(N);
        facteur.Zero();
        facteur.Set(N-1, 0, z[0] + z[N-1]);
        facteur.Set(0, N-1, -z[0] - z[N-1]);

        for(int i=0; i<N-1; i=i+1){
            facteur.Set(i, i+1, z[i] + z[i+1]);
            facteur.Set(i+1, i, -z[i] - z[i+1]);
        }

        MaClasse facteur_ = facteur*(1.0/2)*(1.0/dx);
      
        double coeff = eps/6;
        double coeff_ = eps/2;
        MaClasse M = D1_ + D3_*coeff + facteur*coeff_;
   
        std::vector <double> b(N, 0.0);

        MaClasse Mat = I + M*(-dt/2);
        Mat.Mlt(xi_1, b);
        
        MaClasse A = I + M*(dt/2);
        A.Solve(b, xi_2); // Résolution du système linéaire (I + dt/2M)X_n+1 = (I-dt/2*M)X_n


        for (int k=0; k<x_space.size(); k++){
            results[i+1][k] = xi_1[k];
            xi_0[k] = xi_1[k];
            xi_1[k] = xi_2[k];
        }    
    }


    // Ecriture des résultats dans un fichier .txt pour pour les afficher avec Gnuplot ( une colonne = une solution)
    ofstream fichier("results.txt", ios::out | ios::trunc);
    if (fichier)
    {
        for(int i=0; i< N; i++){
            for (int j=0; j<n_time+1; j++){
                if (j == 0)
                    fichier << x_space[i] << "\t";
                if (j==1)
                    fichier <<results[j-1][i]<< "\t";
                if ((j+1)%500==0)
                    fichier << results[j][i] << "\t";
            }
            fichier << endl;

        }
        fichier.close(); 
    }

    else 
        std::cerr << "Impossible d'ouvrir le fichier"<< endl;    
        
    return 0;

}