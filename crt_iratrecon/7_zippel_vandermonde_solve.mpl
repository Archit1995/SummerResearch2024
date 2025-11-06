Zippel_Vandermonde_solver:=proc(y::list,terms::integer,roots_::list,lambda_,p::integer)# Correct so far. 
    local M,fin_coeff,q,q_lambda_inv,V_inv_b,i,j:
    lprint("In Zippel_Vandermonde_solver"):
    # lprint("y=",y):
    # lprint("terms=",terms):
    # lprint("roots_=",roots_):
    # lprint("lambda_=",lambda_):
    M:=lambda_ mod p:
    # lprint("M=",M):
    # lprint("roots=",roots_):
    fin_coeff:=Vector(terms,0):
    for i from 1 to terms do
        q:=quo(M,Z-roots_[i],Z):
        # lprint("q=",q):
        q_lambda_inv:= 1/ Eval(q,Z=roots_[i]) mod p:
        # lprint("q_lambda_inv=",q_lambda_inv):
        qq:=q*q_lambda_inv mod p;
        # lprint("q*q_lambda_inv mod p=",qq):
        V_inv_b:=0:
        for j from 1 to terms do
            # lprint("j=",j):
            # lprint("coeff(q,Z,j-1)=",coeff(q,Z,j-1)):
            # lprint("y[j]=",y[j]):
            # lprint("coeff(q,Z,j-1)*y[j]=",coeff(q,Z,j-1)*y[j] mod p):
            V_inv_b:=V_inv_b+coeff(q,Z,j-1)*y[j] mod p:
        end do:
        # lprint("V_inv_b=",V_inv_b):
        fin_coeff[i]:=V_inv_b*q_lambda_inv mod p:
    end do:
    # lprint("fin_coeff=",fin_coeff):
    return convert(fin_coeff,list):
end proc: