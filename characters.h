#include "slint.h"

#include <complex>
#include <cmath>
#include <vector>
#include <fftw3.h>

using namespace std;

#include "mpfi_fft.h"

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
    complex<double> sum(long end);
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
 
    mpfi_c_t * zeta_powers_odd_mpfi;
    mpfi_c_t * zeta_powers_even_mpfi;

    mpfi_c_t twopii_on_qodd;
    mpfi_c_t twopii_on_qeven;

    void DFTsum(std::complex<double> * out, std::complex<double> * in);
    void DFTsum(mpfi_c_t * out, mpfi_c_t * in);

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

        k = 0;
        phi_q = euler_phi(q);
        dft_translation_built = false;

        if(q_odd > 1) {
            primes = new std::vector<long>();
            exponents = new std::vector<int>();
            generators = new std::vector<long>();

            factors(q_odd, primes, exponents);
            k = primes->size();
            phi_q_odd = euler_phi(q_odd);
            
            PHI = new long[k];
            A = new long*[q_odd];
            for(int n = 0; n < q_odd; n++) {
                A[n] = new long[k];
            }
            zeta_powers_odd = new std::complex<double>[phi_q_odd];
            zeta_powers_odd_mpfi = new mpfi_c_t[phi_q_odd];

            for(long j = 0; j < k; j++) {
                long x = pow(primes->at(j), exponents->at(j));
                long g = primitive_root(x);
                generators->push_back(g);
                long phi = (pow(primes->at(j), exponents->at(j) - 1) * (primes->at(j) - 1));
                PHI[j] = phi_q_odd/phi;
                long a = 1;
                for(long l = 0; l < phi; l++) {
                    for(long m = a; m < q_odd; m += x) {
                        A[m][j] = l;
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

            mpfi_c_init(twopii_on_qodd);
            mpfi_c_zero(twopii_on_qodd);
            mpfi_const_pi(twopii_on_qodd[0].im);
            mpfi_mul_ui(twopii_on_qodd[0].im, twopii_on_qodd[0].im, 2ul);
            mpfi_div_ui(twopii_on_qodd[0].im, twopii_on_qodd[0].im, phi_q_odd);

            for(unsigned long n = 0; n < phi_q_odd; n++) {
                zeta_powers_odd[n] = e(n/(double)phi_q_odd);
                mpfi_c_init(zeta_powers_odd_mpfi[n]);
                mpfi_c_mul_ui(zeta_powers_odd_mpfi[n], twopii_on_qodd, n);
                mpfi_c_exp(zeta_powers_odd_mpfi[n], zeta_powers_odd_mpfi[n]);
            }
        } // end of initialization of everything having to do
          // with q_odd
        
        if(q_even > 4) {
            B = new long[q_even];
            zeta_powers_even = new complex<double>[q_even/4];
            zeta_powers_even_mpfi = new mpfi_c_t[q_even/4];

            mpfi_c_init(twopii_on_qeven);
            mpfi_c_zero(twopii_on_qeven);
            mpfi_const_pi(twopii_on_qeven[0].im);
            mpfi_mul_ui(twopii_on_qeven[0].im, twopii_on_qeven[0].im, 2ul);
            mpfi_div_ui(twopii_on_qeven[0].im, twopii_on_qeven[0].im, q_even/4);

            for(unsigned long n = 0; n < q_even/4; n++) {
                zeta_powers_even[n] = e(4*n/double(q_even));
                mpfi_c_init(zeta_powers_even_mpfi[n]);
                mpfi_c_mul_ui(zeta_powers_even_mpfi[n], twopii_on_qeven, n);
                mpfi_c_exp(zeta_powers_even_mpfi[n], zeta_powers_even_mpfi[n]);
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
            for(int n = 0; n < phi_q_odd; n++) {
                mpfi_c_clear(zeta_powers_odd_mpfi[n]);
            }
            delete [] zeta_powers_odd_mpfi;
            mpfi_c_clear(twopii_on_qodd);
            delete [] A;
            delete [] PHI;
            delete primes;
            delete exponents;
            delete generators;
        }
        if(q_even > 4) {
            delete [] B;
            delete [] zeta_powers_even;
            for(int n = 0; n < q_even/4; n++) {
                mpfi_c_clear(zeta_powers_even_mpfi[n]);
            }
            delete [] zeta_powers_even_mpfi;
            mpfi_c_clear(twopii_on_qeven);
        }

        if(dft_translation_built) {
            delete [] dft_translation;
            delete [] idft_translation;
            delete [] dft_lengths;
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
        if(m >= q_odd)
            m %= q_odd;
        if(n >= q_odd);
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

    void chi(mpfi_c_t out, long m, long n) {
        complex<double> odd_part = 1.0;
        if(q_even > 1) {
            if(m % 2 == 0 || n % 2 == 0) {
                mpfi_c_zero(out);
                return;
            }
            else if(q_even == 2) {
                mpfi_c_set_ui(out, 1, 0);
            }
            else if(q_even == 4) {
                if(m % 4 == 3 && n % 4 == 3) {
                    mpfi_set_si(out[0].re, -1);
                    mpfi_set_si(out[0].im, 0);
                }
                else {
                    mpfi_c_set_ui(out, 1, 0);
                }
            }
            else {
                    mpfi_c_set(out, zeta_powers_even_mpfi[chi_even_exponent(m % q_even, n % q_even)]);
            }
        }
        else mpfi_c_set_ui(out, 1, 0);

        if(m >= q_odd)
            m %= q_odd;
        if(n >= q_odd);
            n %= q_odd;
        if(q_odd > 1) {
            if(A[m][0] == -1 || A[n][0] == -1)
                mpfi_c_zero(out);
            else
                mpfi_c_mul(out, out, zeta_powers_odd_mpfi[chi_odd_exponent(m, n)]);
        }
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

    void DFTsum_direct (mpfi_c_t * out, unsigned long * in) {
        mpfi_c_t x, y;
        mpfi_c_init(x);
        mpfi_c_init(y);

        for(int n = 0; n < q; n++) {
            mpfi_c_zero(out[n]);
            if(!is_coprime_to_q(n)) continue;
            for(int k = 0; k < q; k++) {
                if(!is_coprime_to_q(k)) continue;
                chi(y, n, k);
                //cout << endl;
                mpfi_c_mul_ui(x, y, in[k]);
                mpfi_c_add(out[n], out[n], x);
            }
        }

        mpfi_c_clear(x);
        mpfi_c_clear(y);
    }

    void all_sums(complex<double> * out, long end) {
        complex<double> * dft_in = new complex<double>[q]();
        for(long k = 0; k <= end; k++) {
            dft_in[k] = 1.0;
        }
        DFTsum(out, dft_in);
        delete [] dft_in;
    }

    void all_sums(mpfi_c_t * out, long end) {
        mpfi_c_t * dft_in = new mpfi_c_t[q];
        for(long k = 0; k < q; k++) {
            mpfi_c_init(dft_in[k]);
            if(k <= end) {
                mpfi_c_set_ui(dft_in[k], 1ul, 0ul);
            }
            else {
                mpfi_c_zero(dft_in[k]);
            }
        }
        DFTsum(out, dft_in);
        for(int k = 0; k < q; k++) {
            mpfi_c_clear(dft_in[k]);
        }
        delete [] dft_in;
    }



    int * dft_translation;
    int * idft_translation;
    int * dft_lengths;
    int dft_dimension;
    bool dft_translation_built;

    void build_dft_translation(bool force=false) {
        if(dft_translation_built && !force) return;
        if(dft_translation_built && force) {
            delete [] dft_translation;
            delete [] idft_translation;
            delete [] dft_lengths;
        }

        int even_dimension;
        if(q_even <= 2)      even_dimension = 0;
        else if(q_even == 4) even_dimension = 1;
        else                 even_dimension = 2;

        dft_dimension = k + even_dimension;
        dft_lengths = new int[dft_dimension];
        
        if(even_dimension >= 1)
            dft_lengths[0] = 2;
        if(even_dimension == 2)
            dft_lengths[1] = q_even/4;

        for(int n = 0; n < k; n++) {
            long p = primes->at(n);
            long e = exponents->at(n);
            dft_lengths[n+even_dimension] = (p - 1)*pow(p, e-1);
        }

        dft_translation = new int[q];
        idft_translation = new int[q];
        int dlog;
        //          m == g[j]**A[m][k] mod p[j]**e[j]
        for(int n = 0; n < q; n++) {
            if((q_even > 1 && n % 2 == 0) || (q_odd > 1 && A[n % q_odd][0] == -1)) {
                dft_translation[n] = -1;
                idft_translation[n] = -1;
            }
            else {
                int dft_index = 0;
                int idft_index = 0;
                if(even_dimension >= 1) {
                    idft_index += (n % 4 == 3);
                    dft_index += (n % 4 == 3);
                    if(dft_dimension > 1) {
                        idft_index *= dft_lengths[1];
                        dft_index *= dft_lengths[1];
                    }
                }
                if(even_dimension == 2) {
                    dlog = B[n % q_even];
                    idft_index += dlog;
                    if(dlog > 0)
                        dft_index += (q_even/4 - B[n % q_even]);
                    if(dft_dimension > 2) {
                        idft_index *= dft_lengths[2];
                        dft_index *= dft_lengths[2];
                    }
                }
                for(int j = 0; j < k; j++) {
                    dlog = A[n % q_odd][j];
                    idft_index += dlog;
                    if(dlog > 0)
                        dft_index += (dft_lengths[j + even_dimension] - A[n % q_odd][j]);
                    if(j < k-1) {
                        idft_index *= dft_lengths[j + even_dimension + 1];
                        dft_index *= dft_lengths[j + even_dimension + 1];
                    }
                }
                //cout << dft_index << " " << idft_index << endl;
                dft_translation[n] = dft_index;
                idft_translation[n] = idft_index;
            }
        }

        dft_translation_built = true;
    }
};

void DirichletGroup::DFTsum (complex<double> * out, complex<double> * in) {
    if(q < 4) {
        // just don't want to deal with problems with 2 right now.
        DFTsum_direct(out, in);
        return;
    }

    build_dft_translation();

    complex<double> *a, *X;
    a = new complex<double>[phi_q];
    X = new complex<double>[phi_q];


    for(int n = 0; n < q; n++) {
        if(dft_translation[n] == -1) continue;
        else a[dft_translation[n]] = in[n];
    }

    fftw_plan plan = fftw_plan_dft(dft_dimension, dft_lengths,
                                        (fftw_complex *)a,
                                        (fftw_complex *)X,
                                        FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);
    for(int n = 0; n < q; n++) {
        if(idft_translation[n] == -1) out[n] = 0.0;
        else out[n] = X[idft_translation[n]];
    }

    fftw_destroy_plan(plan);

}

void DirichletGroup::DFTsum (mpfi_c_t * out, mpfi_c_t * in) {
    //if(q < 4) {
    //    // just don't want to deal with problems with 2 right now.
    //    DFTsum_direct(out, in);
    //    return;
    //}

    build_dft_translation();

    mpfi_c_t *a;//, *X;
    a = new mpfi_c_t[phi_q];
    //X = new mpfi_c_t[phi_q];

    for(int n = 0; n < phi_q; n++) {
        mpfi_c_init(a[n]);
        //mpfi_c_init(X[n]);
    }


    for(int n = 0; n < q; n++) {
        // This might be inefficient.
        // Since the mpfi_c array is really just an array of pointers,
        // we might be able to just copy the pointers without copying
        // the data. The fft wouldn't be operating on continuous memory, then,
        // but it still might be better than this copy.
        if(dft_translation[n] == -1) continue;
        else mpfi_c_set(a[dft_translation[n]], in[n]);
    }

    unsigned long * lengths = new unsigned long[dft_dimension];
    for(int n = 0; n < dft_dimension; n++) lengths[n] = dft_lengths[n];
    ndft(a, phi_q, dft_dimension, lengths);
    //fftw_plan plan = fftw_plan_dft(dft_dimension, dft_lengths,
    //                                    (fftw_complex *)a,
    //                                    (fftw_complex *)X,
    //                                    FFTW_BACKWARD, FFTW_ESTIMATE);

    for(int n = 0; n < q; n++) {
        if(idft_translation[n] == -1) mpfi_c_zero(out[n]);
        else mpfi_c_set(out[n], a[idft_translation[n]]);
    }

    for(int n = 0; n < phi_q; n++) {
        mpfi_c_clear(a[n]);
    }

    delete [] a;

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

complex<double> DirichletCharacter::sum(long end) {
    std::complex<double> S(0,0);

    for(long n = 0; n <= end; n++) {
        S = S + value(n);
    }

    return S;
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

