#pragma once
#define LONG long long int
#include<vector>

extern "C" void polprint64s(LONG* coeffs, int degree, LONG p);

class Polynomial {
    public:
        std::vector<LONG> coefficients;
        std::vector<int> exponents;
        int degree;
        LONG p;
    public:
        Polynomial(std::vector<LONG> c, std::vector<int> e, int deg, LONG p_):coefficients(c),exponents(e),degree(deg),p(p_) {}

        Polynomial(std::vector<LONG>c, int deg,LONG p_):coefficients(c),degree(deg),p(p_){
            exponents.resize(deg+1);
            for(int i=0;i<=deg;i++){
                exponents[i]=i;
            }
        }
        Polynomial():degree(0),p(0){}
        // void print(std::vector<LONG> coefficients, int degree, LONG p){
        void print(){
            polprint64s(const_cast<LONG*>(coefficients.data()), degree, p);
        }
    };



    

