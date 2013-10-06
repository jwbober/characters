#include "slint.h"

#include <complex>
#include <cmath>
#include <vector>
#include <fftw3.h>

using namespace std;

inline std::complex<double> e(double z) {
    std::complex<double>twopii(0,2*M_PI);
    return std::exp(twopii * z);
}

class DirichletGroup;


class DirichletCharacter {
public:
    long m;         // the label
    DirichletGroup * parent;

    DirichletCharacter() {}
    DirichletCharacter(DirichletGroup * parent_, long m_);

    ~DirichletCharacter() {}
    long exponent(long n);
    complex<double> max(long* index);
    std::complex<double> value(long n);
    bool is_primitive();
    bool is_primitive_at_two();
    bool is_primitive_at_odd_part();
    bool is_even();

    complex<double> gauss_sum();
};


class DirichletGroup {
    //
    //
public:
    long q;         // the modulus
    long q_odd;     // the odd part of the modulus
    long q_even;

    int k;          // the number of odd prime factors of q
 
    std::vector<long> * primes;
    std::vector<int> * exponents;
    std::vector<long> * generators;

    long ** A;      // exponent vectors:
                    //      for each m coprime to q we store an array
                    //      with the property that
                    //
                    //          m == g[j]**A[m][k] mod p[j]**e[j]
                    //
                    //      where g[j] is a primitive root mod p[j]**e[j],
                    //      and p[j] is the j-th prime factor of q.
                    //      (We don't actually store g[k], p[j], or e[j] anywhere.)
    
    long * B;

    long * PHI;     // PHI[j] = phi(q/phi(p[j]**e[j])). This will make it easier
                    // to compute the characters.

    long phi_q_odd;
    long phi_q;     // phi(q)


    std::complex<double> * zeta_powers_odd;     // zeta_powers[n] == e(n/phi(q))
    std::complex<double> * zeta_powers_even;     // zeta_powers[n] == e(n/phi(q))
 
    void DFTsum(std::complex<double> * out, std::complex<double> * in);

    bool is_coprime_to_q(long n) {
        if(q_even > 1 && n % 2 == 0)
            return false;
        else if(q_odd == 1)
            return true;
        else
            return A[n % q_odd][0] != -1;
    }

    DirichletGroup(long q_) : q(q_) {
        q_even = 1;
        q_odd = q;
        while(q_odd % 2 == 0) {
            q_odd /= 2;
            q_even *= 2;
        }

        if(q_odd > 1) {
            primes = new std::vector<long>();
            exponents = new std::vector<int>();
            generators = new std::vector<long>();

            factors(q_odd, primes, exponents);
            k = primes->size();
            phi_q_odd = euler_phi(q_odd);
            phi_q = euler_phi(q);
            
            PHI = new long[k];
            A = new long*[q_odd];
            for(int n = 0; n < q_odd; n++) {
                A[n] = new long[k];
            }
            zeta_powers_odd = new std::complex<double>[phi_q_odd];

            for(long j = 0; j < k; j++) {
                long x = pow(primes->at(j), exponents->at(j));
                long g = primitive_root(x);
                generators->push_back(g);
                long phi = (pow(primes->at(j), exponents->at(j) - 1) * (primes->at(j) - 1));
                PHI[j] = phi_q_odd/phi;
                long a = 1;
                //cout << exponents[k] << endl;
                for(long l = 0; l < phi; l++) {
                    for(long m = a; m < q_odd; m += x) {
                        A[m][j] = l;
                        //std::cout << m << std::endl;
                    }
                    a = (a * g) % x;
                }
            }
        
            //
            // for each m, 0 <= m < q, if (m,q) > 1, store
            // this as a flag in A[m][0]
            //
            for(long m = 0; m < q_odd; m++) {
                if(GCD(m,q_odd) > 1) {
                    A[m][0] = -1;
                }
            }

            for(long n = 0; n < phi_q_odd; n++) {
                zeta_powers_odd[n] = e(n/(double)phi_q_odd);
            }
        } // end of initialization of everything having to do
          // with q_odd
        
        if(q_even > 4) {
            B = new long[q_even];
            zeta_powers_even = new complex<double>[q_even/4];
            for(long n = 0; n < q_even/4; n++) {
                zeta_powers_even[n] = e(4*n/double(q_even));
            }
            long pow_five = 1;
            for(long e = 0; e < q_even/4; e++) {
                B[pow_five] = e;
                B[pow_five - 1] = 1;
                B[q_even - pow_five] = e;
                B[q_even - pow_five - 1] = -1;
                pow_five *= 5;
                pow_five %= q_even;
            }
        }

    }

    ~DirichletGroup() {
        if(q_odd > 1) {
            delete [] zeta_powers_odd;
            for(int n = 0; n < q_odd; n++) {
                delete [] A[n];
            }
            delete [] A;
            delete [] PHI;
            delete primes;
            delete exponents;
            delete generators;
        }
        if(q_even > 4) {
            delete [] B;
            delete [] zeta_powers_even;
        }
    }

    long chi_odd_exponent(long m, long n) {
        long x = 0;
        for(int j = 0; j < k; j++) {
            x += A[m][j]*A[n][j]*PHI[j];
            //x = x % phi_q_odd;
            if(x > 4294967296)
                x = x % phi_q_odd;
        }
        if(x >= phi_q_odd) 
            x = x % phi_q_odd;

        return x;
    }
    
    long chi_even_exponent(long m, long n) {
        long x = B[m]*B[n];
        if(B[m-1] == -1 && B[n-1] == -1)
            x += q_even/8;
        return x % (q_even/4);
    }

    std::complex<double> chi(long m, long n) {
        complex<double> even_part = 1.0;
        complex<double> odd_part = 1.0;
        if(q_even > 1) {
            if(m % 2 == 0 || n % 2 == 0) {
                return 0;
            }
            else if(q_even == 2) even_part = 1.0;
            else if(q_even == 4) {
                if(m % 4 == 3 && n % 4 == 3) even_part = -1.0;
                else even_part = 1.0;
            }
            else {
                even_part = zeta_powers_even[chi_even_exponent(m % q_even, n % q_even)];
            }
        }
        if(m > q_odd)
            m %= q_odd;
        if(n > q_odd);
            n %= q_odd;
        if(q_odd > 1) {
            if(A[m][0] == -1 || A[n][0] == -1)
                return 0;
            else
                odd_part = zeta_powers_odd[chi_odd_exponent(m, n)];
            if(q_even == 1)
                return odd_part;
        }
        return odd_part * even_part;
    }

    DirichletCharacter character(long m) {
        return DirichletCharacter(this, m);
    }

    void DFTsum_direct (complex<double> * out, complex<double> * in) {
        //
        // Set out[n] to the sum
        //
        //     sum_{k=1}^{q-1} in[k] chi(n,k)
        //
        // (computed naively, for testing purposes, or because
        //  something better hasn't been implemented yet.)

        
        for(int n = 0; n < q; n++) {
            complex<double> S = 0.0;
            for(int k = 0; k < q; k++) {
                S += in[k] * chi(n,k);
            }
            out[n] = S;
        }
    }
};

void DirichletGroup::DFTsum (complex<double> * out, complex<double> * in) {
    // The easiest case is when q is an odd prime power.
    if(q_even == 1 && k == 1) {
        // the ordering of the character labels is not very FFT
        // friendly, so most of the code here is going to be
        // just to get things in the right order
        
        fftw_complex *a, *X;
        a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(phi_q));
        X = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(phi_q));
        
        fftw_plan plan = fftw_plan_dft_1d(phi_q, a, X,
                                          FFTW_BACKWARD, FFTW_ESTIMATE);

        long g = generators->at(0);
        long z = 1;

        // translate in one direction
        for(int k = 0; k < phi_q; k++) {
            //a[k] = in[z];
            a[k][0] = in[z].real();
            a[k][1] = in[z].imag();
            z = z * g;
            z = z % q;
        }

        fftw_execute(plan);

        // translate back
        z = 1;
        for(int k = 0; k < phi_q; k++) {
            out[z] = complex<double>(X[k][0], X[k][1]);
            z = z * g;
            z = z % q;
        }

        delete [] a;
        delete [] X;
    }
    else if(q_even == 1) {
        complex<double> *a, *X;
        a = new complex<double>[phi_q];
        X = new complex<double>[phi_q];

        int * phi = new int[k];
        for(int n = 0; n < k; n++) {
            long p = primes->at(n);
            long e = exponents->at(n);
            phi[n] = (p - 1)*pow(p, e-1);
        }

        int * translation = new int[q];
        //          m == g[j]**A[m][k] mod p[j]**e[j]
        for(int n = 0; n < q; n++) {
            if(A[n][0] == -1) {
                translation[n] = -1;
            }
            else {
                int index = 0;
                for(int j = 0; j < k; j++) {
                    index += A[n][j];
                    if(j < k-1) index *= phi[j+1];
                }
                translation[n] = index;
            }
        }

        for(int n = 0; n < q; n++) {
            if(translation[n] == -1) continue;
            else a[translation[n]] = in[n];
        }
        
        fftw_plan plan = fftw_plan_dft(k, phi, (fftw_complex *)a,
                                               (fftw_complex *)X,
                                       FFTW_BACKWARD, FFTW_ESTIMATE);
        
        fftw_execute(plan);

        for(int n = 0; n < q; n++) {
            if(translation[n] == -1) continue;
            else out[n] = X[translation[n]];
        }

        delete [] a;
        delete [] X;
    }
}




DirichletCharacter::DirichletCharacter(DirichletGroup * parent_, long m_) : parent(parent_), m(m_) {
        //
        // This doesn't actually do anything. Computing characters
        // using DirichletCharacter is not going to be any faster
        // than using DirichletGroup::chi(), at least for now.
        //
}

complex<double> DirichletCharacter::gauss_sum() {
    complex<double> S = 0;
    long q = parent->q;
    complex<double> z = e(1.0/q);
    complex<double> x = z;
    for(long n = 1; n < q; n++) {
        S = S + value(n) * x;
        x = x * z;
    }
    return S;
}

complex<double> DirichletCharacter::max(long * index) {
    // return the max of the partial sums of the character,
    // and set *index to the spot where the max occurs,
    // unless index = 0
    //
    // note that *index is always <= (q-1)/2
    
    //if(parent->A[m][0] == -1) {
    //    if(index != 0){
    //        *index = -1;
    //    }
    //    return 0;
    //}

    std::complex<double> S(0,0);
    double absmax = 0.0;
    complex<double> current_max = 0.0;
    long max_location;

    for(long n = 0; n <= (parent->q-1)/2; n++) {
        S = S + value(n);
        if(abs(S) > absmax) {
            absmax = abs(S);
            current_max = S;
            max_location = n;
        }
    }

    if(index != 0) {
        *index = max_location;
    }
    return current_max;
}

long DirichletCharacter::exponent(long n) {
    long x = 0;
    for(int j = 0; j < parent->k; j++) {
        x += parent->A[m][j]*parent->A[n][j]*parent->PHI[j];
        if(x > 4294967296)
            x = x % parent->phi_q;
    }
    if(x >= parent->phi_q) {
        x = x % parent->phi_q;
    }
    return x;
}

std::complex<double> DirichletCharacter::value(long n) {
    return parent->chi(m,n);
    //if(parent->A[m][0] == -1 || parent->A[n][0] == -1)
    //    return 0;
    //return parent->zeta_powers[exponent(n)];
}

bool DirichletCharacter::is_primitive() {
    // return whether or not this character is primitive

    return is_primitive_at_two() && is_primitive_at_odd_part();
}

bool DirichletCharacter::is_primitive_at_odd_part() {
    // this is computed one prime at a time, and the
    // character will be primitive if and only if it
    // is primitive at every prime.
    if(parent->q_odd == 1) return true;
    else {
        long n = m % parent->q_odd;
        if(parent->A[n][0] == -1) return false;
        for(int j = 0; j < parent->k; j++) {
            if(parent->A[n][j] % parent->primes->at(j) == 0)
                return false;
        }
        return true;
    }
}

bool DirichletCharacter::is_primitive_at_two() {
    long q_even = parent->q_even;
    long * B = parent->B;
    long n = m % q_even;
    if(q_even == 1) return true;
    else if(q_even == 2) return false;
    else if(q_even == 4) return n == 3;
    else {
        n = n % 8;
        if(n == 3) return true;
        return n == 5;
    }
}

bool DirichletCharacter::is_even() {
    // return whether or not this character is even
    //
    // We just figure out if the character is even by evaluating it at q-1.
    // The evaluation isn't going to be exact, but since the number is just
    // +-1 we can just check which of these it is closest to.
    //
    
    return abs(value(parent->q - 1) - 1.0) < .5;
}

