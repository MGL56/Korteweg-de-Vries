
#include <vector>
#include "maclasse.h"
#include <iostream>
using namespace std;


void MaClasse::Reallocate(int size){

    n = size - 2;
    try{

        if (size < 6)
            std::cout << "Taille N < 6!";
        else{
            d = new double [n];
            e = new double [n-1];
            f = new double [n-2];
            c = new double [n-1];
            b = new double [n-2];
        }
    }
    catch(const std::exception& k){
        d = new double [4];
        e = new double [3];
        f = new double [2];
        c = new double [3];
        b = new double [2];
        std::cout << k.what() << "\n";

    }
}

void MaClasse::Zero(){

    for(int i=0; i<n; i=i+1)
        d[i] = 0;
    

    for(int i=0; i<n-1; i=i+1){
        e[i] = 0;
        c[i] = 0;
    }

    for(int i=0; i<n-2; i=i+1){
        f[i] = 0;
        b[i] = 0;
    }
    
    for (int i=0; i<6; i=i+1){
        row[0] = 0;
        col[0] = 0;
    }

    a[0][0] = 0;
    a[1][0] = 0;
    a[1][1] = 0;
    a[0][1] = 0;

}



void MaClasse::Mlt(std::vector<double> x, std::vector<double> &y){

    int N = n + 2;

    y[0] = x[0]*d[0] + x[1]*e[0] + x[2]*f[0] + col[0]*x[N-2] + col[1]*x[N-1];
    y[1] = x[0]*c[0] + x[1]*d[1] + x[2]*e[1] + x[3]*f[1] + col[2]*x[N-1];

    y[N-2] = row[0]*x[0] + row[3]*x[n-2] + row[4]*x[n-1] + a[0][0]*x[N-2] + a[0][1]*x[N-1];
    y[N-1] = row[1]*x[0] + row[2]*x[1] + row[5]*x[n-3] + a[1][0]*x[N-2] + a[1][1]*x[N-1];

    y[n-2] = b[n-4]*x[n-4] + c[n-3]*x[n-3] + d[n-2]*x[n-2] + e[n-2]*x[n-1] + col[3]*x[N-2];
    y[n-1] = b[n-3]*x[n-3] + c[n-2]*x[n-2] + d[n-1]*x[n-1] + col[4]*x[N-2] + col[5]*x[N-1];

    for(int i=2; i<n-2; i=i+1)
        y[i] = b[i-2]*x[i-2] + c[i-1]*x[i-1] + d[i]*x[i] + e[i]*x[i+1] + f[i]*x[i+2];

}


void MaClasse::FactoriseBand(){

    /**
    Factorisation LU de la matrice bande

    **/

    for (int i=0; i<n-1; i=i+1){
        c[i] = c[i]/d[i];
        d[i+1] -= c[i]*e[i];
        if (i < n-2){
            e[i+1] -= c[i]*f[i];
            b[i] = b[i]/d[i];
            c[i+1] -= b[i]*e[i];
            d[i+2] -= b[i]*f[i];
        }
    }
}

std::vector<double> MaClasse::SolveBand(std::vector<double> x){

    /**
     Résolution d'un système At=x en supposant A factorisée en LU
     <=> Ly = x
         Ut = y

    **/

    std::vector<double> t(x.size(), 0.0), y(x.size(), 0.0);
 

    // Résolution du système Ly=x par descente
    y[0] = x[0];
    y[1] = x[1] - c[0]*x[0];
    for (int i=2; i<x.size(); i++)
         y[i] = x[i] - b[i-2]*x[i-2] - c[i-1]*x[i-1];
    
    // Résolution du système Ut = y par remonté
    t[n-1] = y[n-1]/d[n-1];
    t[n-2] = 1/d[n-2]*(y[n-2] - e[n-2]*t[n-1]);
    for (int i=x.size()-3; i>=0; i--)
        t[i] = 1/d[i]*(y[i] - e[i]*t[i+1] - f[i]*t[i+2]);

    return t;

}

void MaClasse::Factorise(){

    /**
    A = (A0,0 A0,1
         A1,0 A1,1)

    AX = B <=> A0,0*x0 + A0,1*x1 = b0
               A1,0*x0 + A1,1*x1 = b1 

           <=> x0 = A0,0-1(b0 - A0,1*x1)
               (A1,1 - A1,0*A0,0^(-1)*A0,1)x1 = b1 - A1,0*A0,0^(-1)*b0

    **/


    // On factorise la bande en LU
    this->FactoriseBand();

    std::vector<double> b1(n, 0.0), b2(n, 0.0), x2(n, 0.0), x1(n, 0.0);

    b1[0] = col[0];
    b1[n-2] = col[3];
    b1[n-1] = col[4];

    b2[0] = col[1];
    b2[1] = col[2];
    b2[n-1] = col[5];


    x1 = this->SolveBand(b1);
    x2 = this->SolveBand(b2);



    schur[0][0] = a[0][0] - (row[0]*x1[0] + row[3]*x1[n-2] + row[4]*x1[n-1]);
    schur[1][0] = a[1][0] - (row[1]*x1[0] + row[2]*x1[1] + row[5]*x1[n-1]);

    schur[0][1] = a[0][1] - (row[0]*x2[0] + row[3]*x2[n-2] + row[4]*x2[n-1]);
    schur[1][1] = a[1][1] - (row[0]*x2[0] + row[2]*x2[1] + row[5]*x2[n-1]);

}

void MaClasse::Solve(std::vector<double> b, std::vector<double> &x){

    this->Factorise();
    std::vector<double> x_1 = this->SolveBand(b);

    int N = n + 2;
    std::vector<double> result(n, 0.0);
    std::vector <double> y1(2, 0.0), x1(2, 0.0), x2(n, 0.0);

    // Calcul de x1
    y1[0] = b[N-2] - (row[0]*x_1[0] + row[3]*x_1[n-2] + row[4]*x_1[n-1]);
    y1[1] = b[N-1] - (row[1]*x_1[0] + row[2]*x_1[1] + row[5]*x_1[n-1]);
    
        
    double det_schur = schur[0][0]*schur[1][1] - schur[1][0]*schur[0][1];
    //std::cout<<"det_schur" <<det_schur<<endl;

    x1[0] = 1/det_schur*(schur[1][1]*y1[0] - schur[0][1]*y1[1]);
    x1[1] = 1/det_schur*(-schur[1][0]*y1[0] + schur[0][0]*y1[1]);

    // Calcul de x0
    result[0] = x1[0]*col[0] + x1[1]*col[1];
    result[1] = col[2]*x1[1];
    result[n-2] = col[3]*x1[0];
    result[n-1] = col[4]*x1[0] + col[5]*x1[1];

    x2 = this->SolveBand(result);

    for(int i=0; i<x2.size(); i++)
        x[i] = x_1[i] - x2[i];

    x[N-2] = x1[0];
    x[N-1] = x1[1];
    //std::cout<<endl;
    //std::cout << "x[N-2]" << x[N-1] <<endl;
}




void MaClasse::Set(int i, int j, double value){


    while (i > n+1)
        i -= n+2;
    while (i < 0)
        i += n+2;

    while (j > n+1)
        j -= n+2;
    while (j < 0)
        j += n+2;
    
    if ((i < n) && (j < n)){
        if (i == j)
            d[i] = value;
        else if (i == j+1)
            c[j] = value;
        else if (i == j-1)
            e[i] = value;
        else if (i == j-2)
            f[i] = value;
        else if (i == j+2)
            b[j] = value;
        else 
            return;
    }
    else if ((i > n-1) && (j > n-1))
        a[i-n][j-n] = value;
    
    else if ((i<n) && (j>n-1)){
        if ((i==0) && (j>n-1))
            col[j - n] = value;
        else if ((i==1) && (j==n+1))
            col[2] = value;
        else if ((i==n-2) && (j==n))
            col[3] = value;
        else if ((i==n-1) && (j>n-1))
            col[j -(n-4)] = value;
        else
            return;
    }
    else if ((i>n-1) && (j<n)){
        if ((j==0) && (i>n-1))
            row[i - n] = value;
        else if ((j==1) && (i==n+1))
            row[2] = value;
        else if ((j==n-2) && (i==n))
            row[3] = value;
        else if ((j==n-1) && (i>n-1))
            row[i -(n-4)] = value;
        else 
            return;

    }

}


MaClasse MaClasse::operator*(double coeff){

    MaClasse A;
    A.Reallocate(n+2);
    A.Zero();

    for (int i=0; i<n+2; i=i+1){
        for (int j=0; j<n+2; j=j+1)
            A.Set(i, j, this->operator()(i, j)*coeff);
    }

    return A;

}

MaClasse MaClasse::plus(MaClasse A) {

    MaClasse C;
    C.Reallocate(n+2);
    C.Zero();

    double value;
    for (int i=0; i<n+2; i=i+1){
        for (int j=0; j<n+2; j=j+1){
            value = A(i, j) + this->operator()(i, j);
            C.Set(i, j, value);
        }
    }


    return C;
}


MaClasse operator+(MaClasse A, MaClasse B){

    return A.plus(B);

}

double MaClasse::operator()(int i, int j) const{

        
    while (i > n+1)
        i -= n+2;
    while (i < 0)
        i += n+2;
    while (j > n+1)
        j -= n+2;
    while (j < 0)
        j += n+2;

    if ((i < n) && (j < n)){
        if (i==j)
            return d[i];
        else if (i==j+1)
            return c[j];
        else if (i==j-1)
            return e[i];
        else if (i==j+2)
            return b[j];
        else if (i==j-2)
            return f[i];
        else 
            return 0.0;
    }
    else if  ((i > n-1) && (j > n-1))
        return a[i-n][j-n];

    else if ((i<n) && (j>n-1)){
        if ((i==0) && (j>n-1))
            return col[j - n];
        else if ((i==1) && (j==n+1))
            return col[2];
        else if ((i==n-2) && (j==n))
            return col[3];
        else if ((i==n-1) && (j>n-1))
            return col[j -(n-4)];
        else
            return 0.0;
    }
    else if ((i>n-1) && (j<n)){
        if ((j==0) && (i>n-1))
            return row[i - n];
        else if ((j==1) && (i==n+1))
            return row[2];
        else if ((j==n-2) && (i==n))
            return row[3];
        else if ((j==n-1) && (i>n-1))
            return row[i -(n-4)];
        else
            return 0.0;
    }
    else
        return 0.0;
}

void MaClasse::Affichage(){

    for (int i=0; i<n+2; i=i+1){
        for (int j=0; j<n+2; j=j+1)
            std::cout << float(this->operator()(i, j)) << "\t";
        std::cout << endl;
    }
    
} 
