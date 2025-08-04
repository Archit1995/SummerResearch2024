with(LinearAlgebra):
Small_sys:={x1+y1*x2+y2^2*x3-1,y1^3*x1+x2+y2*x3-2,x1-(y1^2-y2)*x2+y2^3*x3-7};
Vars := { seq(   x||i, i=1..nops(Small_sys) )}:
vars := indets(Small_sys) minus Vars:
for i from 1 to nops(Small_sys) do
    print("Small_sys[", i, "] = ", Small_sys[i]); 
    print("______________________________________");
end do:
# Sort seems to be lexicographic thus x10 is before x2,..,x9
print("Vars = ", Vars);
print("vars = ", vars);

A := GenerateMatrix(Small_sys,Vars,augmented=true):
# A := GenerateMatrix(Sys,Vars):
print("A = ", A);
print("A[1] = ", A[1]);
#gg := convert(WW,listlist):
pars := convert(vars,list);
y:=[1,2,3];
soln:=solve(Small_sys,Vars);
simplify(soln[1]);
p:=23;
[seq(eval(soln[i],{y1=11,y2=3}),i=1..3)]mod p;