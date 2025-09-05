#ifndef POLYALG8_H
#define POLYALG8_H
#endif 
#define LONG long long int
#define ULONG unsigned long long int

#define UINT64 unsigned long long
#include "int128g.c"
ULONG seed;
ULONG mult;
/******************************************************************************************/
/*  Zp utilities                                                                          */
/******************************************************************************************/
LONG rand64s(LONG p);;
LONG mul64s(LONG a, LONG b, LONG p);; 
/* c^(-1); mod p assuming 0 < c < p < 2^63 */
LONG modinv64s( LONG c, LONG p );;  
/* a^n mod p assuming 0 <= a < p < 2^63 */
LONG powmod64s( LONG a, LONG n, LONG p );;
/* a^n mod p assuming 0 <= a < p < 2^63 */
LONG powmodP64s( LONG a, LONG n, LONG p, recint P );;

/******************************************************************************************/
/* Array routines                                                                         */
/******************************************************************************************/

LONG * array(LONG n);;
int * arrayint(LONG n);  
void veccopy64s( LONG *u, int n, LONG *v );
void vecfill64s( LONG x, LONG *A, int n );
void vecprint64s( LONG *A, int n );
void vecprint32s( int *A, int n );
/******************************************************************************************/
/* Polynomial routines                                                                    */
/******************************************************************************************/
void polcopy64s( LONG *A, int d, LONG *B );
/* print an array [a0,a1,...,ad] in form ad*x^d+...+a1*x+a0 */
void polprint64s( LONG *A, int d, LONG p ); 
/* check for coefficients are on [0,p); */
void polcheck64s( LONG *A, int d, LONG p, char *x );
// c = a + b mod p
int poladd64s(LONG *a, LONG *b, LONG *c, int da, int db, LONG p);
// a += b mod p
int poladdto64s(LONG *a, LONG *b, int da, int db, LONG p); 
    // c = a-b mod p
int polsub64s(LONG *a, LONG *b, LONG *c, int da, int db, LONG p); 
// a -= b mod p
int polsubfrom64s(LONG *a, LONG *b, int da, int db, LONG p);

// compute A = A - (ax+b); B efficiently
int polsubmul( LONG *A, LONG *B, LONG a, LONG b, int dA, int dB, LONG p );
/* compute gcd(A,B); and put gcd in A and return it's degree */
int polsubmulP( LONG *A, LONG *B, LONG a, LONG b, int dA, int dB, LONG p, recint P ); 

int poldiff64s( LONG *f, int d, LONG *fp, LONG p ); 
// Q = f/(x-alpha);  return the remainder.  Q must be of size d+1
LONG poldiv164s( LONG *f, int d, LONG alpha, LONG *Q, LONG p, recint P );
// M = (x+alpha); M
void polmul164s( LONG *M, int d, LONG alpha, LONG p, recint P ); 
// w.y mod p
LONG dotprod64s( LONG *w, LONG *y, int n, LONG p ); 

/* compute C(x); = A(x);^2 mod p and return deg(C); */
/* we allow C to overwrite A i.e. polsqr64s(A,A,d,p); */
int polsqr64s( LONG * A, LONG * C, int d, LONG p );


/* compute C(x); = A(x); * B(x); mod p and return deg(C); */
/* we allow C to overwrite either A or B i.e. polmul64s(A,B,A,da,db,p); */
int polmul64s( LONG * A, LONG * B, LONG * C, int da, int db, LONG p);


/* compute C(x); = A(x); * B(x); mod p and return deg(C); */
/* we allow C to overwrite either A or B i.e. polmul64s(A,B,A,da,db,p); */
/* Use mulrec64 */
int polmulrec64s( LONG * A, LONG * B, LONG * C, int da, int db, LONG p );

//  Polynomial fma (fused multiply add); C += A*B          
//  Unlike polmul64s, C must be distinct from A and B
int polfma64s( LONG * A, LONG * B, LONG * C, int da, int db, int dc, LONG p);

/* divide A by B and put the remainder and quotient in A */
/* return the degree of the remainder                    */
int poldiv64s( LONG * A, LONG * B, int da, int db, LONG p );

/* divide A by B and put the remainder and quotient in A */
/* return the degree of the remainder                    */
int poldivrec64s( LONG * A, LONG * B, int da, int db, LONG p );

int mulmod64s( LONG *A, LONG *B, LONG *M, int da, int db, int dm,
    LONG *W, LONG *C, int dc, int flag, LONG p );


void polscamul64s( LONG x, LONG *A, int d, LONG p ); 

/* make polynomial in A monic */
void monic64s( LONG *A, int d, LONG p );

/* compute gcd(A,B); and put gcd in A and return it's degree */
/* Both A and B are destroyed */
int polgcd64s( LONG * A, LONG * B, int da, int db, LONG p );


void polgcdext64s( LONG *A, LONG *B, int da, int db,
    LONG *G, LONG *S, LONG *T, int *dG, int *dS, int *dT,
    //LONG *s1, *s2, *t1, *t2, int *ds1, int *ds2, int *dt1, int *dt2,
    LONG *W, LONG p );


/********************************************************************************/


int monic(LONG *A, int m, int *dA, LONG *M, int d, LONG *G, LONG *S, LONG *T, LONG *W, LONG p );

void printmatrix( LONG *A, int m, int n );

void matrixvecmul64s(LONG *A, LONG *u, int d, LONG *v, LONG p );

// A is an m+1 by d matrix, B is n+1 by d matrix, stored in row major order
// A (and B); encodes a polynomial in R[x] where R=Zp[z]/M(z);
// A = sum( sum( A[i*d+j] z^j, j=0..dA[i] ); x^i, i=0..m );, i.e. m = deg(A,x);
// Compute A = monic gcd(A,B); inplace
// LONG *S, LONG *T, LONG *W, // I'm allocating these here now
int alggcd64s(
    LONG *A, int m, int *dA,LONG *B, int n, int *dB,LONG *M, int d, LONG *G,LONG p );


    

// Apply phi to A[i] and B[i], compute their GCD mod M(z); and invert phi 
int alggcdphi64s( 
    LONG *A, int m, // A is an (m+1); x d matrix
    LONG *B, int n, // B is an (n+1); x d matrix
    LONG *M, int d, // M is a monic polynomial of degree d
    LONG *phi, LONG *phiinv, // d x d matrices
    LONG *G, // for zero divisor
    LONG p );

void transpose( LONG *A, int n ); 


// Consider R = Zp[x,y]/<x^2-2,y^2-3> where p=101
// M = z^4-10z^2+1 the minimial polynomial
// row(0,A); = x = 46z + 51z^3 = [0,46,0,51]
// row(1,A); = y = 56z + 50z^3 = [0,56,0,50]
// Let z = x+y, G = GB([x^2-2,y^2-3,z-x-y],plex(x,y,z);); = [M,x(z);,y(z);].
// Old basis is [1,y,x,xy].  New basis is [1,z,z^2,z^3]
// col(0,B); = 1 = [1, 0,0, 0] 
// col(1,B); = y = [0,46,0,51]
// col(2,B); = x = [0,56,0,50]
// col(3,B); = xy = [48,0,51,0] = x y mod M
void phimapping(
    LONG *M, int dm,
    LONG *A, int n, int *D, // row(i,A); is a generator for X[i]
    LONG *B, // Change of basis matrix
    LONG *W, // Working storage
    LONG p );



// S = 1/A mod B
int polmodinv64s( LONG *A, LONG *M, int da, int dm,LONG *G, LONG *S, LONG *W, LONG p );
/* C(x); := A(x);^n mod B(x); mod p;  0<=deg(A);<deg(B); and R must be of size 2*db-1 */
/* If A(x); is not reduced mod B(x); then we first compute C(x); := A(x); mod B(x);   */
int polpowmod64s( LONG * A, LONG n, LONG * B, int da, int db, LONG *C, LONG *R, LONG p );
// Input f in Zp[x] of degree d > 0, a known product of d linear factors.
// Output roots of f in R.
// The input array f is destroyed.
// W is a scratch array of size at least 3*d
void polsplit64s( LONG *f, int d, LONG *R, LONG *W, LONG p );


int polroots64s( LONG * f, int d, LONG * R, LONG *W, LONG p );

// Input sequence a = [a1,a2,a3,...,aN]
// Output polynomial Lambda(x); is written to L
// Uses the half extended Euclidean algorithm

int BerlekampMassey64s( LONG *a, int N, LONG *L, LONG *W, LONG p );

/* R is an array of LONG roots, n , initialized with 0sis the number of roots
L is a buffer that stores ouput, initialized with 0s
W is a buffer that stores , initialized with 0s
p is a prime */
void polLambda64s( LONG *R, int n, LONG *L, LONG *W, LONG p );