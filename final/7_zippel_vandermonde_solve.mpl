Zippel_Vandermonde_solver:=proc(y::list,terms::integer,roots_::list,lambda_,p::integer)# Correct so far. 
    local M,fin_coeff,q,q_lambda_inv,V_inv_b,i,j:
    print("In Zippel_Vandermonde_solver"):
    print("y=",y):
    print("terms=",terms):
    print("roots_=",roots_):
    print("lambda_=",lambda_):
    M:=lambda_ mod p:
    # print("M=",M):
    # print("roots=",roots_):
    fin_coeff:=Vector(terms,0):
    for i from 1 to terms do
        q:=quo(M,Z-roots_[i],Z):
        # print("q=",q):
        q_lambda_inv:= 1/ Eval(q,Z=roots_[i]) mod p:
        # print("q_lambda_inv=",q_lambda_inv):
        V_inv_b:=0:
        for j from 1 to terms do
            V_inv_b:=V_inv_b+coeff(q,Z,j-1)*y[j] mod p:
        end do:
        fin_coeff[i]:=V_inv_b*q_lambda_inv mod p:
    end do:
    # print("fin_coeff=",fin_coeff):
    return convert(fin_coeff,list):
end proc: