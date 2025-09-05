
// This file has my classical O(d^2) arithmetic library for Fp[x]
// It supports 63 bit primes i.e., p < 2^63
// It eliminates the O(d^2) divisions by p using accumulators
// It also uses Roman Pearce's mulrec64 routine for multiplication mod p from the file "int128g.c"
// Copyright Michael Monagan 2000--2022.

// polyalg4.c  inplace polLambda(R,n,L,W,p);
// polyalg5.c  inplace mulmod64s  alggcd64s  phimapping
// polyalg6.c  inplace denserpoly  polfma64s
// polyalg7.c  recden multiplication in R[x] where R=Zp[y]/M(y).  recunpackn
// polyalg8.c  remove recden routines and cleanup

// Compile with gcc -O3 -shared -o polyalg8.so -fPIC polyalg8.c


#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define LONG long long int
#define ULONG unsigned long long int

/******************************************************************************************/
/*  Zp utilities                                                                          */
/******************************************************************************************/

#define UINT64 unsigned long long
// #include "polyalg8.h"
#include "int128g.c"
//typedef struct {
//        UINT64 s;       /* shift */
//        UINT64 v;       /* reciprocal of d */
//        UINT64 d0;      /* divisor shifted up */
//        UINT64 d1;
//} recint;
//
//recint recip1(UINT64  p);
//UINT64 mulrec64(UINT64  a, UINT64  b, recint  v);


ULONG seed;
ULONG mult;

LONG rand64s(LONG p) {
    LONG x,y;
    extern ULONG seed, mult;
    seed = mult*seed;
    x = seed >> 32;
    seed = mult*seed;
    y = seed >> 32;
    x = (x<<31) | y;
    x = x % p;
    return(x);
}


int  min32s(int a, int b) { if( a<b ) return a; else return b; }
int  max32s(int a, int b) { if( a>b ) return a; else return b; }
LONG add64s(LONG a, LONG b, LONG p) { LONG t; t = (a-p)+b; t += (t>>63) & p; return t; }
LONG sub64s(LONG a, LONG b, LONG p) { LONG t; t = a-b; t += (t>>63) & p; return t; }
LONG neg64s(LONG a, LONG p) { return (a==0) ? 0 : p-a; }
LONG mod64s(LONG a, LONG p) { 
    if( a<=-p || a>=p ) a = a%p; 
    if( a<0 ) a += p; 
    return a;
}
LONG mul64s(LONG a, LONG b, LONG p) {
        LONG q, r;
        __asm__ __volatile__(           \
        "       mulq    %%rdx           \n\t" \
        "       divq    %4              \n\t" \
        : "=a"(q), "=d"(r) : "0"(a), "1"(b), "rm"(p));
        return r;
}


        /* z += a1:a0 */
        #define zadd(z,a0,a1) __asm__(\
        "       addq    %4, %0  \n\t" \
        "       adcq    %5, %1  \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]) : "0"(z[0]), "1"(z[1]), "r"(a0), "r"(a1))

        /* z -= a1:a0 */
        #define zsub(z,a0,a1) __asm__(\
        "       subq    %4, %0  \n\t" \
        "       sbbq    %5, %1  \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]) : "0"(z[0]), "1"(z[1]), "r"(a0), "r"(a1))

        /* z = a*b */
        #define zmul(z,a,b) __asm__(\
        "       mulq    %%rdx   \n\t" \
                : "=a"(z[0]), "=d"(z[1]) : "a"(a), "d"(b))

        /* z += a*b */
        #define zfma(z,a,b) do {        \
        unsigned long u,v;              \
        __asm__(                        \
        "       mulq    %%rdx           \n\t" \
        "       addq    %%rax, %0       \n\t" \
        "       adcq    %%rdx, %1       \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]), "=a"(u), "=d"(v) : "0"(z[0]), "1"(z[1]), "a"(a), "d"(b));\
        } while (0)

        /* z -= a*b */
        #define zfms(z,a,b) do {        \
        unsigned long u,v;              \
        __asm__(                        \
        "       mulq    %%rdx           \n\t" \
        "       subq    %%rax, %0       \n\t" \
        "       sbbq    %%rdx, %1       \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]), "=a"(u), "=d"(v) : "0"(z[0]), "1"(z[1]), "a"(a), "d"(b));\
        } while (0)
        /* z[0] = z % p */
        /* z[1] = z / p */
        /* quotient can overflow */
        #define zdiv(z,p) __asm__(\
        "       divq    %4      \n\t" \
                : "=a"(z[1]), "=d"(z[0]) : "a"(z[0]), "d"(z[1]), "r"(p))

        /* z = z % p safe */
        #define zmod(z,p) __asm__(\
        "       divq    %4      \n\t" \
        "       xorq    %0, %0  \n\t" \
                : "=a"(z[1]), "=d"(z[0]) : "a"(z[0]), "d"(z[1] < p ? z[1] : z[1] % p), "r"(p))

        /* z = z << s */
        #define zshl(z,s) __asm__(\
        "       shldq   %%cl, %0, %1    \n\t" \
        "       shlq    %%cl, %0        \n\t" \
                : "=&r"(z[0]), "=&r"(z[1]) : "0"(z[0]), "1"(z[1]), "c"(s))


/* c^(-1) mod p assuming 0 < c < p < 2^63 */
LONG modinv64s( LONG c, LONG p )
{   LONG d,r,q,r1,c1,d1;
    d = p;
    c1 = 1;
    d1 = 0;
    while( d != 0 ) {
        q = c / d;
        r = c - q*d; r1 = c1 - q*d1;
        c = d; c1 = d1;
        d = r; d1 = r1;
    }
    if( c!=1 ) return( 0 );
    if( c1 < 0 ) c1 += p;
    return( c1 );
}


/* a^n mod p assuming 0 <= a < p < 2^63 */
LONG powmod64s( LONG a, LONG n, LONG p )
{   LONG r,s;
    a += (a>>63) & p; // protect from bad input
    if( n==0 ) return 1;
    if( n==1 ) return a;
    for( r=1, s=a; n>0; n /= 2 ) { if( n & 1 ) r = mul64s(r,s,p); s = mul64s(s,s,p); }
    return r;
}

/* a^n mod p assuming 0 <= a < p < 2^63 */
LONG powmodP64s( LONG a, LONG n, LONG p, recint P )
{   LONG r,s;
    a += (a>>63) & p; // protect from bad input
    if( n==0 ) return 1;
    if( n==1 ) return a;
    //for( r=1, s=a; n>0; n /= 2 ) { if( n & 1 ) r = mul64s(r,s,p); s = mul64s(s,s,p); }
    for( r=1, s=a; n>0; n /= 2 ) { if( n & 1 ) r = mulrec64(r,s,P); s = mulrec64(s,s,P); }
    return r;
}

/******************************************************************************************/
/* Array routines                                                                         */
/******************************************************************************************/

LONG * array(LONG n) {//array on heap
    LONG *A;
    A = (LONG *) malloc( n*sizeof(LONG) );
    if( A==0 ) { printf("out of memory\n"); exit(1); }
    return A;
}

int * arrayint(LONG n) {
    int *A;
    A = (int *) malloc( n*sizeof(int) );
    if( A==0 ) { printf("out of memory\n"); exit(1); }
    return A;
}

void veccopy64s( LONG *u, int n, LONG *v ) {
    int i;
    for( i=0; i<n; i++ ) v[i] = u[i];
    return;
}

void vecfill64s( LONG x, LONG *A, int n )
{   int i;
    for( i=0; i<n; i++ ) A[i] = x;
    return;
}

/* print an array in form [a0,a1,...,an-1] */
void vecprint64s( LONG *A, int n )
{   int i;
    printf("[");
    for( i=0; i<n; i++ ) { printf("%lld",A[i]); if( i<n-1 ) printf(", "); }
    printf("]");
    return;
}

void vecprint32s( int *A, int n )
{   int i;
    printf("[");
    for( i=0; i<n; i++ ) { printf("%d",A[i]); if( i<n-1 ) printf(", "); }
    printf("];\n");
    return;
}

/******************************************************************************************/
/* Polynomial routines                                                                    */
/******************************************************************************************/

void polcopy64s( LONG *A, int d, LONG *B )
{   int i;
    for( i=0; i<=d; i++) B[i]=A[i];
    return;
}

/* print an array [a0,a1,...,ad] in form ad*x^d+...+a1*x+a0 */
void polprint64s( LONG *A, int d, LONG p ) {
    int i;
    if( d==-1 ) { printf("0;\n"); return; }
    for( i=d; i>0; i-- ) if( A[i]!=0 ) printf("%lld*x^%d+",A[i],i);
    printf("%lld;\n",A[0]);
    return;
}


/* check for coefficients are on [0,p) */
void polcheck64s( LONG *A, int d, LONG p, char *x ) {
    int i;
    for( i=0; i<=d; i++ )
        if( A[i]<0 || A[i]>=p ) {
            // printf("%s[%d]=%lld\n");
            // printf("%s = "); 
            polprint64s(A,d,p); 
            printf("\n");
            exit(1);
        }
}


int poladd64s(LONG *a, LONG *b, LONG *c, int da, int db, LONG p) {
// c = a + b mod p
   int i,m;
   m = min32s(da,db);
   for( i=0; i<=m; i++ ) c[i] = add64s(a[i],b[i],p);
   if( da==db ) { while( da>=0 && c[da]==0 ) da--; return da; }
   if( da<db ) { for( i=da+1; i<=db; i++ ) c[i] = b[i]; return db; }
   for( i=db+1; i<=da; i++ ) c[i] = a[i]; return da;
}

int poladdto64s(LONG *a, LONG *b, int da, int db, LONG p) {
// a += b mod p
   int i,m;
   m = min32s(da,db);
   for( i=0; i<=m; i++ ) a[i] = add64s(a[i],b[i],p);
   if( da>db ) return da;
   if( da==db ) { while( da>=0 && a[da]==0 ) da--; return da; }
   for( i=da+1; i<=db; i++ ) a[i] = b[i];
   return db;
}


int polsub64s(LONG *a, LONG *b, LONG *c, int da, int db, LONG p) {
// c = a-b mod p
   int i,m;
   m = min32s(da,db);
   for( i=0; i<=m; i++ ) c[i] = sub64s(a[i],b[i],p);
   if( da==db ) { while( da>=0 && c[da]==0 ) da--; return da; }
   if( da>db ) { for ( i=db+1; i<=da; i++ ) c[i] = a[i]; return da; }
   for( i=da+1; i<=db; i++ ) c[i] = neg64s(b[i],p); return db;
}


int polsubfrom64s(LONG *a, LONG *b, int da, int db, LONG p) {
// a -= b mod p
   int i,m;
   m = min32s(da,db);
   for( i=0; i<=m; i++ ) a[i] = sub64s(a[i],b[i],p);
   if( da>db ) return da;
   if( da==db ) { while( da>=0 && a[da]==0 ) da--; return da; }
   for( i=da+1; i<=db; i++ ) a[i] = neg64s(b[i],p);
   return db;
}


int polsubmul( LONG *A, LONG *B, LONG a, LONG b, int dA, int dB, LONG p ) {
   // compute A = A - (ax+b) B efficiently

   LONG t; int i; ULONG z[2];

   if( dB==-1 ) return dA; // B = 0
   z[0] = z[1] = 0LL;

   // if deg(A) <= deg(B) then pad A with zeroes
   while( dA<=dB ) A[++dA] = 0;

   // constant term is special
   t = mul64s(b,B[0],p) ;
   A[0] = sub64s(A[0],t,p);

   for( i=1; i<=dB; i++ ) { zmul(z,a,B[i-1]); zfma(z,b,B[i]); zmod(z,p);
        t = A[i]-z[0]; A[i] = t + ((t>>63)&p); }

   // update leading term from B
   t = mul64s(a,B[dB],p);
   A[dB+1] = sub64s(A[dB+1],t,p);

   // compute and return degree
   while( dA>=0 && (A[dA]==0 || A[dA]==p) ) dA--;
   return dA;
}


/* compute gcd(A,B) and put gcd in A and return it's degree */
int polsubmulP( LONG *A, LONG *B, LONG a, LONG b, int dA, int dB, LONG p, recint P ) {

   // compute A = A - (ax+b) B efficiently

   LONG s,t; int i, d;

   if( dB==-1 ) return dA; // B = 0

   d = dA;

   //polcheck64s(A,dA,p,"A");
   //polcheck64s(B,dB,p,"B");

   // if deg(A) <= deg(B) then pad A with zeroes
   while( dA<=dB ) A[++dA] = 0;

   // constant term is special
   t = mulrec64(b,B[0],P);
   A[0] = sub64s(A[0],t,p);

   //for( i=1; i<=dB; i++ ) { t = mul64s(a,B[i-1],p); t = add64s(t,mul64s(b,B[i],p),p); A[i] = sub64s(A[i],t,p); }
   for( i=1; i<=dB; i++ ) { t = mulrec64(a,B[i-1],P); t = add64s(t,mulrec64(b,B[i],P),p); A[i] = sub64s(A[i],t,p); }

   // update leading term from B
   t = mulrec64(a,B[dB],P);
   A[dB+1] = sub64s(A[dB+1],t,p);

   // compute and return degree
   while( dA>=0 && (A[dA]==0 || A[dA]==p) ) dA--;

   if( dA==d ) printf("polsubmul failure: dAin=%d dAout=%d\n",d,dA);
   return dA;
}


int poldiff64s( LONG *f, int d, LONG *fp, LONG p ) {
    int i; recint P;
    P = recip1(p);
    for( i=1; i<=d; i++ ) fp[i-1] = mulrec64(f[i],(LONG) i,P);
    for( d--; d>=0 && fp[d]==0; d-- );
    return d;
}

LONG poldiv164s( LONG *f, int d, LONG alpha, LONG *Q, LONG p, recint P ) {
    // Q = f/(x-alpha)  return the remainder.  Q must be of size d+1
    int i; LONG r;
    //for( i=0; i<=d; i++ ) Q[i] = f[i];
    //for( i=d; i>0; i-- ) Q[i-1] = add64s(Q[i-1],mulrec64(Q[i],alpha,P),p);
    Q[d] = f[d];
    for( i=d; i>0; i-- ) Q[i-1] = add64s(f[i-1],mulrec64(Q[i],alpha,P),p);
    r = Q[0];
    for( i=0; i<d; i++ ) Q[i] = Q[i+1];
    return r; // remainder
}


void polmul164s( LONG *M, int d, LONG alpha, LONG p, recint P ) {
    // M = (x+alpha) M
    int i;
    for( i=d; i>=0; i-- ) M[i+1] = M[i];
    M[0] = 0;
    for( i=0; i<=d; i++ ) M[i] = add64s(M[i],mulrec64(M[i+1],alpha,P),p);
    return;
}

//printf("u:="); vecprint64s(u,n);
//printf("v:="); vecprint64s(v,n);
    //s = 0;

LONG dotprod64s( LONG *w, LONG *y, int n, LONG p ) 
{   // w.y mod p
    int i; LONG s; ULONG z[2];
    //for( i=0,s=0; i<n; i++ ) s = add64s(s,mulrec64(w[i],y[i],P),p);
    z[0] = z[1] = 0;
    for( i=0; i<n-1; ) {
        zfma(z,w[i],y[i]); i++; 
        zfma(z,w[i],y[i]); i++;
        if( z[1]>=p ) z[1] -= p;
    }
    if( i==n-1 ) zfma(z,w[i],y[i]);    
    zmod(z,p);
    s = z[0];
    return s;
}


/* compute C(x) = A(x)^2 mod p and return deg(C) */
/* we allow C to overwrite A i.e. polsqr64s(A,A,d,p) */
int polsqr64s( LONG * A, LONG * C, int d, LONG p )
{
    int i,k,m,dc; ULONG z[2];
    if( d<0 ) return d;
    for( k=2*d; k>=0; k-- ) {
       m = min32s(k,d);
       i = max32s(0,k-d);
       z[0] = z[1] = 0ll;
       while( i<m-2 ) {
            zfma(z,A[i++],A[m--]);
            if( z[1]>=p ) z[1] -= p;
            zfma(z,A[i++],A[m--]);
       }
       if( i<m ) {
            zfma(z,A[i++],A[m--]);
            if( z[1]>=p ) z[1] -= p;
       }
       zadd(z,z[0],z[1]);
       if( z[1]>=p ) z[1] -= p;
       if( i==m ) zfma(z,A[i],A[i]);
       zmod(z,p);
       C[k] = z[0];
    }
    for( dc = 2*d; dc>=0 && C[dc]==0; dc-- );
    // Why is this loop here? Z_p has no zero-divisors.
    // Because p may not be prime!!
    return( dc );
}


/* compute C(x) = A(x) * B(x) mod p and return deg(C) */
/* we allow C to overwrite either A or B i.e. polmul64s(A,B,A,da,db,p) */
int polmul64s( LONG * A, LONG * B, LONG * C, int da, int db, LONG p)
{
    int i,k,m;
    if( da<0 || db<0 ) return -1;
    int dc = da+db;
if( p<2147483648ll ) { LONG t, p2;
    p2 = p<<32;
    for( k=dc; k>=0; k-- ) {
       i = max32s(0,k-db);
       m = min32s(k,da);
       t = 0ll;
       while( i<m ) {
           t -= A[i]*B[k-i]; i++;
           t -= A[i]*B[k-i]; i++;
           t += (t>>63) & p2;
       }
       if( i==m ) t -= A[i]*B[k-i];
       t = (-t) % p;
       t += (t>>63) & p;
       C[k] = t;
    }
} else { ULONG z[2];
    for( k=dc; k>=0; k-- ) {
       i = max32s(0,k-db);
       m = min32s(k,da);
       z[0] = z[1] = 0ll;
       while( i<m ) {
           zfma(z,A[i],B[k-i]); i++;
           zfma(z,A[i],B[k-i]); i++;
           if( z[1]>=p ) z[1] -= p;
       }
       if( i==m ) zfma(z,A[i],B[k-i]);
       zmod(z,p);
       C[k] = z[0];
    }
}
    for( ; dc>=0 && C[dc]==0; dc-- );
    return( dc );
}


/* compute C(x) = A(x) * B(x) mod p and return deg(C) */
/* we allow C to overwrite either A or B i.e. polmul64s(A,B,A,da,db,p) */
/* Use mulrec64 */
int polmulrec64s( LONG * A, LONG * B, LONG * C, int da, int db, LONG p )
{
    int i,k,m; LONG t;
    recint P = recip1(p);
    if( da<0 || db<0 ) return -1;
    int dc = da+db;
    for( k=dc; k>=0; k-- ) {
       i = max32s(0,k-db);
       m = min32s(k,da);
       for( t=0; i<=m; i++ ) t = add64s(t, mulrec64(A[i],B[k-i],P), p );
       C[k] = t;
    }
    for( ; dc>=0 && C[dc]==0; dc-- );
    return( dc );
}


//  Polynomial fma (fused multiply add) C += A*B          
//  Unlike polmul64s, C must be distinct from A and B
int polfma64s( LONG * A, LONG * B, LONG * C, int da, int db, int dc, LONG p)
{   // polynomial fma (fused multiply add) C += A*B
    int i,k,m; ULONG z[2];
    if( da<0 || db<0 ) return dc; // dc not -1
    for( k=0; k<=da+db; k++ ) {
       i = max32s(0,k-db);
       m = min32s(k,da);
       z[0] = z[1] = 0ll;
       while( i<m ) {
           zfma(z,A[i],B[k-i]); i++;
           zfma(z,A[i],B[k-i]); i++;
           if( z[1]>=p ) z[1] -= p;
       }
       if( i==m ) zfma(z,A[i],B[k-i]);
       zmod(z,p);
       if( k>dc ) C[k] = z[0]; else C[k] = add64s(C[k],z[0],p);
    }
    for( dc=max32s(dc,da+db); dc>=0 && C[dc]==0; dc-- );
    return( dc );
}


/* divide A by B and put the remainder and quotient in A */
/* return the degree of the remainder                    */
int poldiv64s( LONG * A, LONG * B, int da, int db, LONG p )
{
    int dq,dr,k,j,m; LONG t,inv;
    if( db<0 ) { printf("division by zero\n"); exit(1); }
    if( da<db ) return da; 

    if( da==db && B[db]==1 ) {
        // A = [c,b,a] = ax^2+bx+c  B = [D,C,1] = x^2+Cx+D
        // Update A to be [c-aD,b-aC,a].
        t = A[da];
        for( k=0; k<da; k++ )
            if( B[k] ) A[k] = sub64s(A[k],mul64s(t,B[k],p),p);
        for( dr=da-1; dr>=0 && A[dr]==0; dr-- );
        return dr;
    };

    dq = da-db;
    dr = db-1;
    if( B[db]==1 ) inv = 1; else inv = modinv64s(B[db],p);
if( p<2147483648ll ) { LONG p2;
    p2 = p<<32;
    for( k=da; k>=0; k-- ) {
        t = A[k];
        m = min32s(dr,k);
        j = max32s(0,k-dq);
        while( j<m ) {
            t -= B[j]*A[k-j+db]; j++;
            t -= B[j]*A[k-j+db]; j++;
            t += (t>>63) & p2;
        }
        if( j==m ) t -= B[j]*A[k-j+db];
        t = t % p;
        t += (t>>63) & p;
        if( k>=db && inv!=1 ) t = mul64s(t,inv,p);
        A[k] = t;
    }
} else { ULONG z[2];
    for( k=da; k>=0; k-- ) {
        z[0] = z[1] = 0ll;
        m = min32s(dr,k);
        j = max32s(0,k-dq);
        while( j<m ) {
            zfma(z,B[j],A[k-j+db]); j++;
            zfma(z,B[j],A[k-j+db]); j++;
            if( z[1]>=p ) z[1] -= p;
        }
        if( j==m ) zfma(z,B[j],A[k-j+db]);
        zmod(z,p);
        t = A[k] - z[0];
        t += (t>>63) & p;
        if( k>=db && inv!=1 ) t = mul64s(t,inv,p);
        A[k] = t;
    }
}
    while( dr>=0 && A[dr]==0 ) dr--;
    return( dr );
}


/* divide A by B and put the remainder and quotient in A */
/* return the degree of the remainder                    */
int poldivrec64s( LONG * A, LONG * B, int da, int db, LONG p )
{
    int dq,dr,k,j,m; LONG t,inv;
    recint P = recip1(p);
    if( db<0 ) { printf("division by zero\n"); exit(1); }
    if( da<db ) return da; 
    dq = da-db;
    dr = db-1;
    if( B[db]==1 ) inv = 1; else inv = modinv64s(B[db],p);
    for( k=da; k>=0; k-- ) {
        t = A[k];
        m = min32s(dr,k);
        j = max32s(0,k-dq);
        for( t=A[k]; j<=m; j++ )
            t = sub64s(t,mulrec64(B[j],A[k-j+db],P),p);
        if( k>=db && inv!=1 ) t = mulrec64(t,inv,P);
        A[k] = t;
    }
    while( dr>=0 && A[dr]==0 ) dr--;
    return( dr );
}


int mulmod64s( LONG *A, LONG *B, LONG *M, int da, int db, int dm,
     LONG *W, LONG *C, int dc, int flag, LONG p )
// if flag=0 then C = A*B mod M
// if flag=1 then C += A B mod M
// if flag=-1 then C -= A B mod M
// return deg C
{    int dr,i;
     while( da>=0 && A[da]==0 ) da--;
     while( db>=0 && B[db]==0 ) db--;
     dr = polmul64s(A,B,W,da,db,p);
     dr = poldiv64s(W,M,dr,dm,p);
     // printf("mulmod64s: da=%d db=%d dc=%d flag=%d dr=%d\n",da,db,dc,flag,dr);
//printf("Cbef = "); polprint64s(C,dc,p);
//printf("W = "); polprint64s(W,dr,p);
     if( flag==0 ) {
          for(i=0; i<=dr; i++ ) C[i] = W[i];
          for(i=dr+1; i<dm; i++ ) C[i] = 0; // zero this out
          return dr;
     }
     if( flag==1 ) dc = poladd64s(C,W,C,dc,dr,p); else dc = polsub64s(C,W,C,dc,dr,p);
     while( dc>=0 && C[dc]==0 ) dc--;
//printf("Caft = "); polprint64s(C,dc,p);
     return dc;
}


void polscamul64s( LONG x, LONG *A, int d, LONG p ) {
    int i;
    if( x==1 ) return;
    if( x==-1 ) for( i=0; i<=d; i++ ) A[i] = neg64s(A[i],p);
    else for( i=0; i<=d; i++ ) A[i] = mul64s(x,A[i],p);
    return;
}


/* make polynomial in A monic */
void monic64s( LONG *A, int d, LONG p ) {
    int i; LONG inv;
    if( d<0 || A[d]==1 ) return;
    inv = modinv64s(A[d],p);
    for( i=0; i<d; i++ ) A[i] = mul64s(inv, A[i], p);
    A[d] = 1;
    return;
}


/* compute gcd(A,B) and put gcd in A and return it's degree */
/* Both A and B are destroyed */
int polgcd64s( LONG * A, LONG * B, int da, int db, LONG p ) {
    int dr; LONG *C, *D, *R, u, a, b;
    recint P;
    if( db<0 ) { printf("division by zero\n"); exit(1); }
    P = recip1(p);
    C = A; D = B;
    if( da<db ) { R = C; C = D; D = R; dr = da; da = db; db = dr; }
    while( 1 ) {
        if( db>0 && da-db==1 ) { // normal case
            u = modinv64s(D[db],p);
            a = mulrec64(C[da],u,P);
        //    a = mul64s(C[da],u,p);
            b = mulrec64(a,D[db-1],P);
        //    b = mul64s(a,D[db-1],p);
            b = mulrec64(u,sub64s(C[da-1],b,p),P);  // quotient = a x + b
        //    b = mul64s(u,sub64s(C[da-1],b,p),p);  // quotient = a x + b
        //    dr = polsubmulP(C,D,a,b,da,db,p,P);  // C = C - (a x + b) D
            dr = polsubmul(C,D,a,b,da,db,p);  // C = C - (a x + b) D
            if( dr>=db ) printf("failure\n");
        }
        else dr = poldiv64s(C,D,da,db,p);
        if( dr<0 ) { /* D|C so gcd(A,B)=D */
            if( D!=A ) polcopy64s(D,db,A);
            monic64s( A, db, p );
            return db;
        }
        R = C; C = D; D = R; da = db; db = dr;
        //printf("da=%d db=%d\n",da,db);
    }
}


void polgcdext64s( LONG *A, LONG *B, int da, int db,
                  LONG *G, LONG *S, LONG *T, int *dG, int *dS, int *dT,
                  //LONG *s1, *s2, *t1, *t2, int *ds1, int *ds2, int *dt1, int *dt2,
                  LONG *W, LONG p )
{
    // Solve S A + T B = G = monic gcd(A,B) for G,S,T in Zp[x]
    // The arrays A and B are used for the remainder sequence so they are destroyed
    // G,S,T must all be of size max(da+1,db+1)
    // if( S==0 ) W is working storage of size max(da+1,db+1)
    // if( T==0 ) W is working storage of size max(da+1,db+1)
    // if S==0 or T==0 then S (and/or T) are not computed

    int m,dr,ds,dt,dq,ds1,ds2,dt1,dt2; LONG a,b,u;
    LONG *q,*r,*r1,*r2,*s,*s1,*s2,*t,*t1,*t2;

    recint P; P = recip1(p);

    if( da<0 || db<0 ) { printf("inputs must be non-zero\n"); exit(1); }
    m = max32s(da+1,db+1);
    r1 = A; r2 = B;
    if(S) { s1 = S; s2 = W;   s1[0]=1; ds1=0; ds2=-1; }
    if(T) { t1 = T; 
            if(S) t2 = W+m; else t2 = W;
            t2[0]=1; dt2=0; dt1=-1;
    }
    while( 1 ) {
        if( db>0 && da-db==1 ) { // normal case
            u = modinv64s(r2[db],p);
            a = mul64s(r1[da],u,p);
            b = mul64s(a,r2[db-1],p);
            b = mul64s(u,sub64s(r1[da-1],b,p),p);             // quotient = a x + b
            //dr = polsubmul(r1,r2,a,b,da,db,p);                // r1 = r1 - (a x + b) r2
            //if(S) ds = polsubmul(s1,s2,a,b,ds1,ds2,p);        // s1 = s1 - (a x + b) s2
            //if(T) dt = polsubmul(t1,t2,a,b,dt1,dt2,p);        // t1 = t1 - (a x + b) t2
            dr = polsubmulP(r1,r2,a,b,da,db,p,P);              // r1 = r1 - (a x + b) r2
            if(S) ds = polsubmulP(s1,s2,a,b,ds1,ds2,p,P);      // s1 = s1 - (a x + b) s2
            if(T) dt = polsubmulP(t1,t2,a,b,dt1,dt2,p,P);      // t1 = t1 - (a x + b) t2
        }
        else {
            dr = poldiv64s(r1,r2,da,db,p);                 // r1 = [remainder,quotient]
            q  = r1+db; dq = da-db;
            if(S) ds = polmul64s(q,s2,G,dq,ds2,p);
            if(S) ds = polsub64s(s1,G,s1,ds1,ds,p);        // s1 = s1 - q s2
            if(T) dt = polmul64s(q,t2,G,dq,dt2,p);
            if(T) dt = polsub64s(t1,G,t1,dt1,dt,p);        // t1 = t1 - q t2
        }
        if( dr<0 ) { /* D|C so gcd(A,B)=D */
            polcopy64s(r2,db,G);
            if(S) if( s2!=S ) polcopy64s(s2,ds2,S);
            if(T) if( t2!=T ) polcopy64s(t2,dt2,T);
            if( G[db]!=1 ) {
                u = modinv64s(G[db],p);
                polscamul64s(u,G,db,p);
                if(S) polscamul64s(u,S,ds2,p);
                if(T) polscamul64s(u,T,dt2,p);
            }
            dG[0] = db;
            if(S) dS[0] = ds2;
            if(T) dT[0] = dt2;
            return;
        }
        r = r1; r1 = r2; r2 = r;  da = db;   db = dr;
        if(S) { s = s1; s1 = s2; s2 = s; ds1 = ds2; ds2 = ds; }
        if(T) { t = t1; t1 = t2; t2 = t; dt1 = dt2; dt2 = dt; }
    }
}


/********************************************************************************/


int monic( LONG *A, int m, int *dA, LONG *M, int d, 
           LONG *G, LONG *S, LONG *T, LONG *W, LONG p ) {
    int i,j,dr,dG,dS,dT;
    for( i=0; i<=dA[m]; i++ ) T[i] = A[m*d+i];
    for( i=0; i<=d; i++ ) T[d+i] = M[i];
    polgcdext64s(T,T+d,dA[m],d,G,S,0,&dG,&dS,&dT,W,p);
    if( dG>0 ) return dG;
    for( i=0; i<m; i++ ) {
         dr = polmul64s(A+i*d,S,W,dA[i],dS,p);
         dr = poldiv64s(W,M,dr,d,p);
         polcopy64s(W,dr,A+i*d);
         dA[i] = dr;
    }
    dA[m] = 0;
    A[m*d] = 1;
    return 0; // dG = 0
}


void printmatrix( LONG *A, int m, int n ) {
    int i,j;
    for( i=0; i<m; i++ ) {
        printf("[%lld",A[i*n]);
        for( j=1; j<n; j++ ) printf(",%lld",A[i*n+j]);
        printf("]\n");
    }
    return;
}


void matrixvecmul64s(LONG *A, LONG *u, int d, LONG *v, LONG p )
{   int i;
    for( i=0; i<d; i++ ) v[i] = dotprod64s(A+i*d,u,d,p); 
    return;
}


int alggcd64s(
        LONG *A, int m, int *dA,
        LONG *B, int n, int *dB,
        // A is an m+1 by d matrix, B is n+1 by d matrix, stored in row major order
        // A (and B) encodes a polynomial in R[x] where R=Zp[z]/M(z)
        // A = sum( sum( A[i*d+j] z^j, j=0..dA[i] ) x^i, i=0..m ), i.e. m = deg(A,x)
        // Compute A = monic gcd(A,B) inplace
        LONG *M, int d, LONG *G, 
        // LONG *S, LONG *T, LONG *W, // I'm allocating these here now
        LONG p )
{
    LONG *S,*T,*W,*TT;
    int i,j,k,dg,dr,*D;

    // S,T,W are temporary arrays for the extended E.A.
    S = array(d+1); T = array(2*(d+1)); W = array(3*(d+1));

    if( m<n ) { k=m; m=n; n=k; TT=A; A=B; B=TT; D=dA; dA=dB; dB=D; }

    dg = monic(A,m,dA,M,d,G,S,T,W,p);
    if( dg>0 ) { free(S); free(T); free(W); return -dg; }
    while( n >= 0 ) {
        dg = monic(B,n,dB,M,d,G,S,T,W,p);
        if( dg>0 ) { free(S); free(T); free(W); return -dg; }
        //dA[m] = -1;
        for( j=m-n,k=0; k<n; k++,j++ ) {
            dA[j] = polsub64s(A+j*d,B+k*d,A+j*d,dA[j],dB[k],p);
        }
        for( m--; m>=0 && dA[m]<0; m-- ); // m = deg(A,x)
        while( m >= n ) {
            for( j=m-n,k=0; k<n; k++,j++ ) {
                if( dB[k]<0 ) continue;
                dr = polmul64s(A+m*d,B+k*d,W,dA[m],dB[k],p);
                dr = poldiv64s(W,M,dr,d,p);
                dA[j] = polsub64s(A+j*d,W,A+j*d,dA[j],dr,p);
            }
            for( m--; m>=0 && dA[m]<0; m-- );
        }
        k=m; m=n; n=k; TT=A; A=B; B=TT; D=dA; dA=dB; dB=D;
    }
    // n < 0 implies B = 0 so A is the GCD
    for( i=0; i<=m; i++ ) { // copy A into B
         dB[i] = dA[i];
         polcopy64s(A+i*d,dA[i],B+i*d);
    }
    free(S); free(T); free(W);
    return m;
}


int alggcdphi64s( 
        LONG *A, int m, // A is an (m+1) x d matrix
        LONG *B, int n, // B is an (n+1) x d matrix
        LONG *M, int d, // M is a monic polynomial of degree d
        LONG *phi, LONG *phiinv, // d x d matrices
        LONG *G, // for zero divisor
        LONG p )
{       int i,j,*dA,*dB,dg;
        LONG *u,*v;
    // Apply phi to A[i] and B[i], compute their GCD mod M(z) and invert phi 
//printf("phi="); printmatrix(phi,d,d);
//printf("p=%lld\n",p);
//printf("A=");
//printmatrix(A,m+1,d);
    dA = arrayint(m+1);
    v = array(d);
    for( i=0,u=A; i<=m; i++,u+=d ) { // u = row(i,A)
        matrixvecmul64s( phi, u, d, v, p ); // v = phi . u^T
        veccopy64s(v,d,u); // for( j=0; j<d; j++ ) u[j] = v[j]; 
        for( j=d-1; j>=0 && u[j]==0; j-- ); // j = deg(u)
        dA[i] = j;
    }
//printf("phiA=");
//printmatrix(A,m+1,d);
    dB = arrayint(n+1);
    for( i=0,u=B; i<=n; i++,u+=d ) {
        matrixvecmul64s( phi, u, d, v, p );
        veccopy64s(v,d,u); // for( j=0; j<d; j++ ) u[j] = v[j]; 
        for( j=d-1; j>=0 && u[j]==0; j-- );
        dB[i] = j;
    }
    dg = alggcd64s(A,m,dA,B,n,dB,M,d,G,p);
    if( dg<0 ) { free(dB); free(dA); free(v); return dg; }
    for( i=0,u=A; i<=dg; i++,u+=d ) {
        for( j=dA[i]+1; j<d; j++ ) u[j] = 0;
        matrixvecmul64s( phiinv, u, d, v, p ); 
        veccopy64s(v,d,u); // for( j=0; j<d; j++ ) u[j] = v[j];
    }
    free(dA); free(dB); free(v);
    return dg;
}

void transpose( LONG *A, int n ) {
    int i,j;
    LONG t;
    for( i=0; i<n; i++ )
        for( j=i+1; j<n; j++ )
            { t = A[i*n+j]; A[i*n+j] = A[j*n+i]; A[j*n+i] = t; }
    return;
}

void phimapping(
    LONG *M, int dm,
    LONG *A, int n, int *D, // row(i,A) is a generator for X[i]
    LONG *B, // Change of basis matrix
    LONG *W, // Working storage
    LONG p )
// Consider R = Zp[x,y]/<x^2-2,y^2-3> where p=101
// M = z^4-10z^2+1 the minimial polynomial
// row(0,A) = x = 46z + 51z^3 = [0,46,0,51]
// row(1,A) = y = 56z + 50z^3 = [0,56,0,50]
// Let z = x+y, G = GB([x^2-2,y^2-3,z-x-y],plex(x,y,z)) = [M,x(z),y(z)].
// Old basis is [1,y,x,xy].  New basis is [1,z,z^2,z^3]
// col(0,B) = 1 = [1, 0,0, 0] 
// col(1,B) = y = [0,46,0,51]
// col(2,B) = x = [0,56,0,50]
// col(3,B) = xy = [48,0,51,0] = x y mod M
{
    int i,j,k,m,da,db,dr;
    clock_t st,pt,tt;
    vecfill64s(0,B,dm*dm); // B = 0
    B[0] = 1;
    m = 1;
    st = clock();
    for( k=n-1; k>=0; k-- ) {
        for( da=dm-1; da>=0 && A[k*dm+da]==0; da-- );
        for( i=1; i<D[k]; i++ ) {
            for( j=0; j<m; j++ ) {
                //printf("k=%d i=%d j=%d\n",k,i,j);
                for( db=dm-1; db>=0 && B[((i-1)*m+j)*dm+db]==0; db-- );
                dr = polmul64s(B+((i-1)*m+j)*dm,A+k*dm,W,db,da,p);
                dr = poldiv64s(W,M,dr,dm,p);
                polcopy64s(W,dr,B+(i*m+j)*dm);
            }
        }
        m = m*D[k];
    }
    pt = clock()-st;
    st = clock();
    transpose(B,dm);
    tt = clock()-st;
    printf("phi time=%.3fs  trans=%.3fs\n",pt/1000000.,tt/1000000.);
    return;
}

// S = 1/A mod B
int polmodinv64s( LONG *A, LONG *M, int da, int dm,
                  LONG *G, LONG *S, LONG *W, LONG p )
{   int dG, dS, dT;
    // W must be able to hold three polynomials of degree dm
    while( da>=0 && A[da]==0 ) da--;
    if( da<0 ) { printf("A is zero\n"); exit(1); }
    while( dm>=0 && M[dm]==0 ) dm--;
    if( da>=dm ) { printf("deg(A) < deg(M) error\n"); exit(1); }
    //printf("dA=%d dM=%d\n",da,dm);
    // copy A and M into W
    polcopy64s(A,da,W); A = W; W += da+1;
    polcopy64s(M,dm,W); M = W; W += dm+1;
    polgcdext64s( A, M, da, dm, G, S, 0, &dG, &dS, &dT, W, p );
    //printf("dG=%d  dS=%d\n",dG,dS);
    if( dG>0 ) return -dG; else return dS;
}


/* C(x) := A(x)^n mod B(x) mod p;  0<=deg(A)<deg(B) and R must be of size 2*db-1 */
/* If A(x) is not reduced mod B(x) then we first compute C(x) := A(x) mod B(x)   */
int polpowmod64s( LONG * A, LONG n, LONG * B, int da, int db, LONG *C, LONG *R, LONG p )
{
    int dc,k,b[63];

    if( n==0 ) { C[0] = 1; return 0; }
    if( da>=db ) da = poldiv64s(A,B,da,db,p);                   // reduce A mod B first
    for( k=0; n>0; k++ ) { b[k]=n&1; n=n/2; }
    polcopy64s(A,da,C);
    dc = da;
    k--;
    while( k>0 ) { k--;
       // Main step: compute C := C^2 mod B in Zp[x]
       //dc = polmul64s(C,C,R,dc,dc,p);                           //printf("deg(R) = %d; R = ",dc); polprint64s(R,dc);
       dc = polsqr64s(C,R,dc,p);                                //printf("deg(R) = %d; R = ",dc); polprint64s(R,dc);
       dc = poldiv64s(R,B,dc,db,p);
       polcopy64s(R,dc,C);                                      //printf("deg(C) = %d; C = ",dc); polprint64s(C,dc);
       if( b[k]==1 ) {                                          //printf(" b[%d]=%d \n", k, b[k] );
           dc = polmul64s(A,C,R,da,dc,p);                       //printf("deg(R) = %d; R = ",dc); polprint64s(R,dc);
           dc = poldiv64s(R,B,dc,db,p);
           polcopy64s(R,dc,C);                                  //printf("deg(C) = %d; C = ",dc); polprint64s(C,dc);
       }
    }
    return dc;
}


// Input f in Zp[x] of degree d > 0, a known product of d linear factors.
// Output roots of f in R.
// The input array f is destroyed.
// W is a scratch array of size at least 3*d
void polsplit64s( LONG *f, int d, LONG *R, LONG *W, LONG p )
{
   int da,dg; LONG alpha, A[2];
   if( d==1 ) { alpha = p-f[0]; R[0] = alpha; return; }
   alpha = rand64s(p); A[1] = 1; A[0] = alpha;
   da = polpowmod64s( A, (p-1)/2, f, 1, d, W, W+d, p );
   if( da==0 ) return polsplit64s(f,d,R,W,p);      // alpha is unlucky, try again
   W[0] = add64s(W[0],1,p);                        // W = (x+alpha)^((p-1)/2) + 1 mod f
   polcopy64s( f, d, W+d );
   dg = polgcd64s( W, W+d, da, d, p );             // g = gcd( W, f ) in W
   if( dg==0 ) return polsplit64s(f,d,R,W,p);      // g = 1 ==> alpha is unlucky, try again
   poldiv64s(f,W,d,dg,p);                          // compute quotient q = f/g destroying f
   polcopy64s(W,dg-1,f);                           // f = [ g mod x^dg followed by q ]
   polsplit64s(f+dg,d-dg,R,W,p);
   f[dg] = 1;
   polsplit64s(f,dg,R+d-dg,W,p);
   return;
}


int polroots64s( LONG * f, int d, LONG * R, LONG *W, LONG p )
{
   int i, da, dg; LONG A[2]; extern ULONG seed,mult;
   clock_t st,et;
printf("roots: deg(f)=%d\n",d);
    // printf("f := "); polprint64s(f,d);
   for( i=0; i<d && f[i]==0; i++ );
   if( i>0 ) { R[0]=0; return( 1 + polroots64s(f+i,d-i,R+1,W,p) ); }
   if( f[d]!=1 ) monic64s(f,d,p);
   A[1] = 1;
   A[0] = 0;
   st = clock();
   da = polpowmod64s( A, p-1, f, 1, d, W, W+d, p );    // W = x^(p-1) mod f
   et = clock();
   printf("Roots: powmod: x^(p-1) mod f = x^%d + ... time = %10lu ms\n", da, (et-st)/1000 );
   //printf("da = %d, a := ",da); polprint64s(W,da);
   if( da==0 && W[0]==1 ) dg = d; // f is all linear factors
   else { W[0] = sub64s(W[0],1,p); dg = polgcd64s( f, W, d, da, p ); }   // f = gcd(f,W-1)
   //printf("g := "); polprint64s(f,dg);
   printf("Roots: def(f)=%d  #roots=%d\n",d,dg);
   if( dg==0 ) return 0;
   seed = 1;
   mult = 6364136223846793003ll;
   st = clock();
   polsplit64s( f, dg, R, W, p );
   et = clock();
   printf("Roots: split time=%10lu ms\n", (et-st)/1000 );
   return dg; // number of roots in R
}


int BerlekampMassey64s( LONG *a, int N, LONG *L, LONG *W, LONG p )
{
    // Input sequence a = [a1,a2,a3,...,aN]
    // Output polynomial Lambda(x) is written to L
    // Uses the half extended Euclidean algorithm
    int i,m,n,dr,dq,dr0,dr1,dv0,dv1,dt;
    LONG *r,*q,*r0,*r1,*v0,*v1,*t,u,A,b;
    //recint P;
    while( N>0 && a[N-1]==0 ) N--; // ignore leading zeroes
    n = N/2;
    N = 2*n;
    if( N==0 ) return -1;
    m = N-1;
    // W is space for r0 = x^N and r1 of degree m and v0 and v1 of degree at most n
    r0 = W; r1 = r0+N+1; v0 = r1+N; v1 = v0+n+1;
    vecfill64s(0,r0,N); r0[N] = 1; dr0 = N;             // r0 = x^(2*n)
    for(i=0; i<N; i++) r1[i] = a[m-i];
    for(dr1=m; dr1>=0 && r1[dr1]==0; dr1--);            // r1 = sum(a[m-i]*x^i,i=0..m)
    if( dr1==-1 ) return -1;
    dv0 = -1;                                           // v0 = 0
    v1[0] = 1; dv1 = 0;                                 // v1 = 1
    //P = recip1(p);
    while( n <= dr1 ) {
        if( dr1>0 && dr0-dr1==1 ) { // normal case
            u = modinv64s(r1[dr1],p);
            A = mul64s(r0[dr0],u,p);
            b = mul64s(A,r1[dr1-1],p);
            b = mul64s(u,sub64s(r0[dr0-1],b,p),p);             // quotient q = A x + b
            //dr = polsubmulP(r0,r1,A,b,dr0,dr1,p,P);            // r0 = r0 - (A x + b) r1
            dr = polsubmul(r0,r1,A,b,dr0,dr1,p);            // r0 = r0 - (A x + b) r1
            //dt = polsubmulP(v0,v1,A,b,dv0,dv1,p,P);            // v0 = v0 - (A x + b) v1
            dt = polsubmul(v0,v1,A,b,dv0,dv1,p);            // v0 = v0 - (A x + b) v1
        } else {
           dr = poldiv64s(r0,r1,dr0,dr1,p);
           q = r0+dr1; dq = dr0-dr1;                           // q = quo(r0,r1)
           dt = polmul64s(q,v1,L,dq,dv1,p);
           dt = polsub64s(v0,L,v0,dv0,dt,p);
        }
        r = r0; r0 = r1; r1 = r; dr0 = dr1; dr1 = dr;         // r0,r1 = r1,rem(r0,r1)
        t = v0; v0 = v1; v1 = t; dv0 = dv1; dv1 = dt;         // v0,v1 = v1,v0 - q*v1
        //printf("r0 = "); polprint64s(r0,dr0);
        //printf("r1 = "); polprint64s(r1,dr1);
        //printf("v0 = "); polprint64s(v0,dv0);
        //printf("v1 = "); polprint64s(v1,dv1);
    }
    if( dv1>=0 ) {
        polcopy64s(v1,dv1,L);
        monic64s(L,dv1,p);
    }
    return dv1;
}

// recursive algorithm
void polLambda64s( LONG *R, int n, LONG *L, LONG *W, LONG p ) {
// Compute L = proc_{i=0}^{n-1} (x-R[i])
// L must be of length n+1 and W must be of length n
// You can use R for W as in polLambda64s( R, n, L, R, p ) which will destroy R
    printf("In polLambda64s \n");
    int i,m,d;
    //termination condition
    if( n==0 ) { printf("in n = 0 \n"); L[0] = 1; printf("L[0]=%d \n",1); return; }
    if( n==1 ) { printf("in n = 1 \n");  L[0] = neg64s(R[0],p);printf("L[0]=%lld \n",L[0]); L[1] = 1;
        printf("L[1]=%d \n",1); return; }
    if( n==2 ) {
        printf("in n =2 \n");
        L[0] = mul64s(R[0],R[1],p);
        L[1] = neg64s(add64s(R[0],R[1],p),p); 
        L[2] = 1;
        printf("L[0] = %lld \n",L[0]);
        printf("L[1] = %lld \n",L[1]);
        printf("L[2] = %d \n",1);
        return; 
    }
    m = n/2; d = n-m; 
    printf("m:= %d \n",m);
    printf("d:= %d \n",d);
    polLambda64s( R, m, L, W, p );           // L = [ x, x, 1, -, - ]  if n=4
    polLambda64s( R+m, d, L+m, W, p );       // L = [ x, x, y, y, 1 ]
    for( i=0; i<n; i++ ) W[i] = L[i];        // W = [ x, x, y, y ]
    polmul64s( W, W+m, L, m-1, d-1, p );     // L = [ a, b, c, -, - ]
    L[n-1] = 0; L[n] = 1;                    // L = [ a, b, c, 0, 1 ]
    poladd64s( W, L+d, L+d, m-1, m-1, p );   //   + [ -, -, y, y, - ] 
    poladd64s( W+m, L+m, L+m, d-1, d-1, p ); //   + [ -, -, x, x, - ]
    return;
}


/*  Test polLambda  */

// int main() {
//     LONG *R, *W, *L, p;
//     int i,n;
//     p = 5;
//     //  p = (p<<55)+1;
//     printf("p := %lld;\n",p);
//     // n = 10;
//     n=2;
//     R = array(n); for( i=0; i<n; i++ ) R[i] = 2*i+1;//odd numbers
//     printf("R := "); vecprint64s(R,n); printf("\n");
//     L = array(n+1);
//     W = array(n);
//     printf("W := %llu \n",W[0]);
//     polLambda64s(R,n,L,W,p);
//     printf("L := "); polprint64s(L,n,p); printf("\n");
//     return 0;
// }


/**** Test polfma64s
int main() {
    LONG *A, *B, *C, p;
    int i,n,d;
    p = 5; p = (p<<55)+1;
    printf("p := %lld;\n",p);
    n = 4;
    A = array(n); for( i=0; i<n; i++ ) A[i] = 2*i+1;
    printf("A := "); polprint64s(A,n-1,p); printf("\n");
    B = array(n); for( i=0; i<n; i++ ) B[i] = 2*i+3;
    printf("B := "); polprint64s(B,n-1,p); printf("\n");
    C = array(2*n); for( i=0; i<2*n; i++ ) C[i] = 0;
    d = polfma64s(A,B,C,n-1,n-1,-1,p);
    d = polfma64s(A,B,C,n-1,n-1,d,p);
    printf("C := "); polprint64s(C,d,p); printf("\n");
    return 0;
}
*****/
