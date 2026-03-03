with(LinearAlgebra):
Constuct_Sys_Blackbox:=proc(Sys,Vars,params) 
    local Lin_BB,L:
    L := GenerateMatrix(Sys,Vars,augmented=true):
    Lin_BB := proc( point_::list(integer), p::prime ) local A,L,T;
        # uses LinearAlgebra:-Modular;
        with(LinearAlgebra[Modular]):
        # L := [[y[1],y[2],1],[y[1]*y[2],-1,1]] mod p;
        L := GenerateMatrix(Sys,Vars,augmented=true) mod p;
        print("L = ",L):
        subs_values := zip((par, pnt) -> par = pnt, params, point_);
        # print("subs_values = ", subs_values):
        L_eval := map(eval, L, subs_values);
        # print("L_eval = ",L_eval):
        num_eqn:=numelems(Vars);
        A := Matrix(num_eqn,num_eqn+1,datatype=integer[8],L_eval);
        print("A = ",A):
        T := traperror( LinearSolve(p,A,1) );
        print(A);
        if T="matrix is singular" then return FAIL fi;
        soln:=convert(A[1..num_eqn,num_eqn+1],list);
        # print("soln = ",soln);
       return soln: # the solution vector x
    end:
end:
# Sys := { y1*x1+y2*x2 - 1, y1*y2*x1-x2 - 1 };
Sys:={x1+y1*x2+y1-3,y2*x1+x2+y1-1}:
#  Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1, 
#                     x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17, 
#                     -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
#                     -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2, 
#                     -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13, 
#                     -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
#                     2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21, 
#                     2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7, 
#                     2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9, 
#                     -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
#                     2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16, 
#                     2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:
p := 2^31-1;
Vars := {seq(   x||i, i=1..nops(Sys) )}:
params := indets(Sys) minus Vars:
params := convert(params,list):
print("params = ", params);
B:=Constuct_Sys_Blackbox(Sys,Vars,params):
y := [3,7,11,13,17];
B(y,p);

