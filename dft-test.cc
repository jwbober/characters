#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

void timesum(long qstart, long qend) {
    complex<double> S = 0.0;
    for(int q = qstart; q < qend; q++) {
        cout << q << endl;
        complex<double> * sums = new complex<double>[q];
        DirichletGroup G(q);
        
        G.all_sums(sums, q/3);
        S += G.chi(2,3);
        delete [] sums;
    }
    cout << S << endl;
}

int main() {
    //timesum( 1 << 14, 1 << 15);
    int qstart = 1;
    double maxerror = 0;
    for(int j = 0; j < 10000; j++) {
        int q = qstart + j;
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

        complex<double> * sums1 = new complex<double>[q];
        complex<double> * sums2 = new complex<double>[q]();

        G.all_sums(sums1, q/3);

        for(int k = 0; k < q; k++) {
            if(G.is_coprime_to_q(k)) {
                DirichletCharacter chi = G.character(k);
                sums2[k] = chi.sum(q/3);
            }
        }
        error = 0;
        for(int k = 0; k < q; k++) {
            error = max(error, abs(sums1[k] - sums2[k]));
        }
        maxerror = max(error, maxerror);
        cout << q << " " << error << " " << maxerror << endl;

        delete [] a;
        delete [] X1;
        delete [] X2;
        delete [] sums1;
        delete [] sums2;
    }
    return 0;
}
