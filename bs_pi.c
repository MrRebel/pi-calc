//#define USE_MINI_GMP

#include "math.h"
#include "stdio.h"
#include "stdint.h"

#ifdef USE_MINI_GMP
#include "gmp-5.1.3/mini-gmp/mini-gmp.c"
#else
#include "gmp.h"
#endif

typedef uint_fast32_t BigInt;
const char* PI_FILE = "./pi.txt";

mpz_t c3over24;
void get_big_sqrt(BigInt, mpz_t);

void bs(BigInt a, BigInt b, mpz_t Pab, mpz_t Qab, mpz_t Tab) {

    mpz_t temp;
    mpz_init(temp);

    if (b == a + 1) { //if they are 1 apart, can calculate exactly
        if (a == 0) { //special case for 0
            mpz_set_ui(Pab, 1);
            mpz_set_ui(Qab, 1);
        } else {

            ////Pab = (6*a-5)*(2*a-1)*(6*a-1)
            BigInt a1, a3;

            a1 = 6*a-5;
            a3 = 6*a-1;

            mpz_set_ui(Pab, 2*a-1);
            mpz_mul_ui(Pab, Pab, a1);
            mpz_mul_ui(Pab, Pab, a3);

            ////Qab = a*a*a*C3_OVER_24
            mpz_set_ui(Qab, a);
            mpz_pow_ui(Qab, Qab, 3);
            mpz_mul(Qab, Qab, c3over24);
        }
        mpz_init_set_ui(Tab, a); // Ta..a+1 = Pa..a+1 * (13591409+54514013*a)
        mpz_mul_ui(Tab, Tab, 545140134);
        mpz_add_ui(Tab, Tab, 13591409);
        mpz_mul(Tab, Tab, Pab);
        if (a&1) { // alternating series ftw7
            mpz_neg(Tab, Tab);
        }
    } else {
        //midpoint of a and b
        BigInt m = (b-a) / 2 + a;

        /// I *could* make new variables buuttt...
        /// I already have Pab, Qab, and Tab!
        /// soooo:
#define Pmb Pab
#define Qam Qab
#define Tam Tab
        mpz_t Pam /*, Qam, Tam*/;
        mpz_t  /*Pmb,*/ Qmb, Tmb;

        mpz_init(Pam);
        // mpz_init(Qam);
        // mpz_init(Tam);
        // mpz_init(Pmb);
        mpz_init(Qmb);
        mpz_init(Tmb);

        /// recursively calculate
        bs(a, m, Pam, Qam, Tam);//use as (a, m)

        bs(m, b, Pmb, Qmb, Tmb);

        mpz_mul(Pab, Pam, Pmb); //multiply

        mpz_mul(Qab, Qam, Qmb);

        mpz_mul(Tam, Tam, Qmb);
        mpz_mul(Tmb, Tmb, Pam);
        mpz_add(Tab, Tam, Tmb);
        

        /// free
        mpz_clear(Pam);
        // mpz_clear(Qam);
        // mpz_clear(Tam);
        // mpz_clear(Pmb);
        mpz_clear(Qmb);
        mpz_clear(Tmb);

#undef Pam
#undef Qam
#undef Tam
    }
    //free the memory!!!
    mpz_clear(temp);
}

void gmp_pi_binary_splitter(BigInt digits, mpz_t pi) {

    mpz_init_set_str(c3over24, "26dd041d878000", 16); // = 10939058860032000

    BigInt N;
    mpz_t P;
    mpz_t Q;
    mpz_t T;
    long double digitsperterm = 14.18164746272547765552552;
        //log10((long double) mpz_get_d(c3over24)/72);

    N = (BigInt) (digits/digitsperterm) + 1;

    mpz_init(P);
    mpz_init(Q);
    mpz_init(T);

    bs(0, N, P, Q, T);

    //Now, P, Q, and T are all set, so we can do stuff (P doesn't really matter for this part)
    mpz_clear(P);

    mpz_t sqrtC;
    mpz_init(sqrtC);

    get_big_sqrt(digits, sqrtC);

    mpz_mul(sqrtC, sqrtC, Q);
    mpz_clear(Q);
    // mpz_mul(sqrtC, sqrtC, base);
    // mpz_mul(sqrtC, sqrtC, base);

    mpz_tdiv_q(pi, sqrtC, T);
    mpz_abs(pi, pi);

    //clear as much memory as possible, then set, clear, and return

    // debug("./T.txt", T);
    // debug("./sqrtC.txt", sqrtC);

    mpz_clear(T);
    mpz_clear(c3over24);//to not cause problems if called again

    mpz_clear(sqrtC);

}

/// TODO: possibly implement Newtons method using a precomputed file of digits for sqrt(10005) to 1M decimal places or so
/// Test against GMP algorithm
void get_big_sqrt(BigInt decimalDigits, mpz_t sqrtC) {
    

    //To get to the final division we have to calculate a huge-ass square root (sqrt(10005)) to the precision "One"
    //where "One" is 10^(# of digits)

// #ifdef USE_MINI_GMP

    mpz_set_ui(sqrtC, 10);
    mpz_pow_ui(sqrtC, sqrtC, decimalDigits * 2);
    mpz_mul_ui(sqrtC, sqrtC, 10005);
    mpz_sqrt(sqrtC, sqrtC);

// #else
    // mp_bitcnt_t binaryPrecision = (mp_bitcnt_t)(3.32192809489 * decimalDigits);
    
    // mpf_t base;
    // mpf_t fsqrtC;

    /// set up the base (maybe make this base 2 later for faster(?) arithmetic TODO: run tests to find out
    // mpf_init2(base, binaryPrecision);
    // mpf_set_ui(base, 10);
    // mpf_pow_ui(base, base, decimalDigits);

    /// get floating point sqrt (much faster!!!)
    // mpf_init2(fsqrtC, binaryPrecision);
    // mpf_sqrt_ui(fsqrtC, 10005);

    /// bring back to an integer value
    // mpf_mul(fsqrtC, fsqrtC, base);

    // mpz_set_f(sqrtC, fsqrtC);
    // mpf_clear(fsqrtC);

// #endif

   mpz_mul_ui(sqrtC, sqrtC, 426880);

}

int main(int argc, char* argv[]) {
    unsigned long long int decimalplaces;
    printf("How many digits of pi?: ");
    while(!scanf("%llu", &decimalplaces))
        printf("Invalid Entry\nHow many digits of pi?: ");

    printf("output goes to %s\n", PI_FILE);

    FILE * pifile;
    pifile = fopen(PI_FILE, "w");
    if (pifile == NULL) {
        fprintf(stderr, "Can't open output file %s!\n", PI_FILE);
        return 1;
    }

    mpz_t pi;
    mpz_init(pi);

    gmp_pi_binary_splitter(decimalplaces, pi);

    mpz_out_str(pifile, 10, pi);

    fclose(pifile);
    mpz_clear(pi);
    return 0;
}
