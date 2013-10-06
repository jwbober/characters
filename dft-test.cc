#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

int main() {
    int qstart = 3;
    double maxerror = 0;
    for(int j = 0; j < 10000; j++) {
        int q = qstart + 2*j;
        DirichletGroup G(q);
        complex<double> * a = new complex<double>[q]();
        complex<double> * X1 = new complex<double>[q]();
        complex<double> * X2 = new complex<double>[q]();

        for(int k = 0; k < q; k++) {
            a[k] = complex<double>(rand()/(double)RAND_MAX, rand()/(double)RAND_MAX);
        }

        G.DFTsum(X1, a);
        G.DFTsum_direct(X2, a);
        
        double error = 0;
        for(int k = 0; k < q; k++) {
            error = max(error, abs(X1[k] - X2[k]));
        }
        maxerror = max(error, maxerror);
        cout << q << " " << error << " " << maxerror << endl;

        delete [] a;
        delete [] X1;
        delete [] X2;
    }
    return 0;
}
