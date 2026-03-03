########################################################
# 0. CONSTRUCT BLACK BOX - For Linear Systems
########################################################
with(LinearAlgebra):
Constuct_Sys_Blackbox:=proc(Sys,Vars,params) 
    local Lin_BB,L:
    L := GenerateMatrix(Sys,Vars,augmented=true):
    Lin_BB := proc( point_::list(integer), p::prime ) 
        local A,T,subs_values,num_eqn,soln,L_eval:
        uses LinearAlgebra:-Modular:
        global counter:
        counter:=counter+1:
        # print("L = ",L):
        subs_values := zip((par, pnt) -> par = pnt, params, point_):
        # print("subs_values = ", subs_values):
        # L_eval := map(eval, L, subs_values):
        L_eval := Eval(L,subs_values) mod p:
        num_eqn:=numelems(Vars):
        A := Matrix(num_eqn,num_eqn+1,datatype=integer[8],L_eval):
        # A := Matrix(num_eqn,num_eqn+1,L_eval):
        # print("Matrix before Gaussian elimination:  ",A):
        T := traperror( LinearSolve(p,A,1) ):
        if T="matrix is singular" then return FAIL fi:
        # print("Matrix after Gaussian elimination: ",A):
        soln:=convert(A[1..num_eqn,num_eqn+1],list):
        # print("soln = ",soln):
       return soln: # the solution vector x
    end:
end:


Construct_Rational_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        local var,num,denom_,a,v:
        global counter:
        counter:=counter+1:       
        var:=vars:
        num:=f:
        denom_:=g:
        a:=num/denom_:
        if Eval(denom_,{seq(var[v]=point_[v],v=1..numelems(point_))}) mod p = 0 then
            # denominator is zero
            error "Denominator is zero":
        else 
            return Eval(a,{seq(var[v]=point_[v],v=1..numelems(point_))}) mod p:
        end if:
    end proc:
    return BB:
end proc:



