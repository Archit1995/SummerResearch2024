with(LinearAlgebra):
with(ArrayTools):

prime_var_map:=table([2=x,3=y,5=z]);
# parser
generate_monomials:=proc(roots_,prime_var_map)
     local ff,l,l2,i;
     monomials:=Vector(numelems(roots_),0);
     for j from 1 to numelems(roots_) do 
          print(roots_[j]);
          ff:=ifactor(roots_[j]);
          print(ff);
          l:=nops(ff);
          for i from 1 to l do
               l2:=nops(op(i,ff)):
               if l2=1 then 
                    ff:=subs(op(i,ff)=prime_var_map[op(1,op(i,ff))],ff);
               else 
                    ff:=subs(op(1,op(i,ff))=prime_var_map[op(1,(op(1,op(i,ff))))],ff);
               fi;
          end do;
          monomials[j]:=ff;
     end do;
     return convert(monomials,list);
end proc;


vars:={x,y,z}:
f:=randpoly(vars,degree=5):
deg:=degree(f):
T:=deg+2:
v:=[seq(eval(f,{x=2^i,y=3^i,z=5^i}),i=0..2*T-1)];
whattype(v):
V:=[seq(v[i..i+(T-1)],i=1..T)]:
H:=Matrix(V);
terms:=Rank(H):
H:=H[1..terms,1..terms]:
S:=-Vector(v[terms+1..terms+terms]):
X:=LinearSolve(H,S):
row:=Size(X)[1]:
Lambda:=z^row:
for i from 1 to row do
    Lambda:=Lambda+X[i]*z^(i-1):
end do:
Lambda;
R:=roots(Lambda):
Roots_:=[seq(r[1],r in R)];
Monomials:=generate_monomials(Roots_,prime_var_map);
