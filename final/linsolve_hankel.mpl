with(LinearAlgebra):
M:=Matrix([[9,12,8],[12,8,9],[8,9,8]]):
b:=Vector([9,8,1]):
p:=17:
c:=LinearSolve(M,b) mod p;
# ch:=sum()
# Roots(c) mod p;