Zippel_Vandermonde_solver:=proc(y::list,terms::integer,roots_::list,lambda_,p::integer)
    local M,fin_coeff,q,q_lambda_inv,V_inv_b,i,j:
    lprint("In Zippel_Vandermonde_solver"):
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
    return convert(fin_coeff,list):
end proc:

