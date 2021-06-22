#include <vector>

class MaClasse{

    friend MaClasse operator+(MaClasse a, MaClasse b);

    private:
        double  *b=NULL, *c=NULL, *d=NULL, *e=NULL, *f=NULL;
        int n;
        double row[6], col[6];
        double a[2][2];
        double schur[2][2];

    public:
        void Zero();
        void Reallocate(int size);
        void Mlt(std::vector<double> x, std::vector<double> &y);
        void Set(int i, int j, double value);
        MaClasse operator*(double);
        double operator() (int, int) const;
        MaClasse plus(MaClasse);
        void Affichage();
        void FactoriseBand();
        std::vector<double> SolveBand(std::vector<double>);
        void Factorise();
        void Solve(std::vector<double> b, std::vector<double> &x);
        


};