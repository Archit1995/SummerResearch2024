#########################################################
# 7. Zippel Transpose Vandermonde Solver
#########################################################

Zippel_Transpose_Vandermonde_solver:=proc(y::list,terms::integer,roots_::list,lambda_::polynom,p::integer)
    description "Solves the Zippel Transpose Vandermonde system to find the coefficients of the polynomial":
# 1. Parameters:
# -----------
# y : list
#     A list of integers representing the evaluations of the univariate images of the numerator or denominator at  powers of 2.
# terms : integer
#     The number of terms (degree + 1) in the polynomial to be reconstructed.
# roots_ : list
#     A list of integers representing the roots of the minimal polynomial obtained from the BMEA algorithm.
# lambda_ : polynom
#     The minimal characteristic polynomial obtained from the BMEA algorithm.
# p : integer
#     A prime number representing the modulus for arithmetic operations.
#
# 2. Returns:
# --------
# list
#     A list of coefficients of the polynomial reconstructed from the Vandermonde system in Z_p.
    local M,fin_coeff,q,q_lambda_inv,V_inv_b,i,j:
    lprint("In Zippel_Transpose_Vandermonde_solver"):
    M:=lambda_ mod p:
    fin_coeff:=Vector(terms,0):
    for i from 1 to terms do
        q:=quo(M,Z-roots_[i],Z):
        q_lambda_inv:= 1/ Eval(q,Z=roots_[i]) mod p:
        V_inv_b:=0:
        for j from 1 to terms do
            V_inv_b:=V_inv_b+coeff(q,Z,j-1)*y[j] mod p:
        end do:
        fin_coeff[i]:=V_inv_b*q_lambda_inv mod p:
    end do:
    # print("final_coeff: ",fin_coeff):
    return convert(fin_coeff,list):
end proc:

