#include <iostream>
#include<vector>
#define LONG long long int
#include "Polynomial.cpp"
extern "C" {
  #include "polyalg8.c"
}

int main() {
    // LONG *roots;
    LONG *W;
    // LONG *coefficients;
    LONG p;
    
    int i,n;
    // p = 2^32-1;
    p=5;
    //  p = (p<<55)+1;
    // printf("p := %lld;\n",p);
    std::cout<<"p := "<<p<<"\n"; 
    // n = 10;
    n=2;
    
    // R = array(n); for( i=0; i<n; i++ ) R[i] = 2*i+1;//odd numbers
    std::vector<LONG>roots(n);
    std::vector<LONG>coefficients(n+1);
    // roots=array(n);
    roots[0]=1;
    roots[1]=1;
    
    // for( i=0; i<n; i++ ) roots[i] =i;

    std::cout<<"roots := \n",
    vecprint64s(const_cast<LONG*>(roots.data()),n); 
    std::cout<<"\n";
    // coefficients = array(n+1);
    W = array(n);
    std::cout<<"coefficients:= "<<coefficients[0]<<"\n";
    // printf("W := %llu \n",W[0]);
    std::cout<<"W:= "<<W[0]<<"\n";
    
    polLambda64s(const_cast<LONG*>(roots.data()),n,const_cast<LONG*>(coefficients.data()),W,p);/* if roots=[1,1] polyLambdas64s
    returns (x-1)(x-1) mod p; 
    */
   Polynomial poly(coefficients,n,p);
    printf("coefficients := ");    
    poly.print();
    



    // polprint64s(const_cast<LONG*>(coefficients.data()),n,p);
    /*Has variable Hardcoded into it.
    Takes array of coefficients,degree and prime p as input */
     printf("\n");
    return 0;
}

