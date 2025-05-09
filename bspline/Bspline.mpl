with(LinearAlgebra):

Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1, x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17, -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
 -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2, -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13, -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
  2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21, 2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7, 2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9, -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
   2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16, 2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:
Vars := { seq(   x||i, i=1..nops(Sys) )}:
vars := indets(Sys) minus Vars:
for i from 1 to nops(Sys) do
    print("Sys[", i, "] = ", Sys[i]); 
    print("______________________________________");
end do:
# Sort seems to be lexicographic thus x10 is before x2,..,x9
print("Vars = ", Vars);
print("vars = ", vars);

A := GenerateMatrix(Sys,Vars,augmented=true):
# A := GenerateMatrix(Sys,Vars):
print("A = ", A);
print("A[1] = ", A[1]);
#gg := convert(WW,listlist):
pars := convert(vars,list);
with(LinearAlgebra):
B := ReducedRowEchelonForm(A):
print("B = ", B);
print("B[1] = ", B[1]);
# x is the last column of B, 22nd entry for each row
x := B[1..21,22]:
print("x = ", x);
n := 21;
for i from 1 to n do factor(x[i]) od:

C := A[1..21,1..21]:
print("C = ", C);
print("C[1] = ", C[1]);
local D := Determinant(C):
#  b contains the rhs of each equation, 22nd entry of each row of A
b := A[1..21,22]:
print("b = ", b):
# print("b[1] = ", b[1]):
print("<C[1..21,1..20]|b> = ", <C[1..21,1..20]|b>);
N := Determinant( <C[1..21,1..20]|b> );
factor(N/D);
for i from 1 to nops(Vars) do 
    print("x",i,"= ",factor(x[i]));
    print("______________________________________");
end do: