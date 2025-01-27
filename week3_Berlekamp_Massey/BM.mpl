

# Maple code for the Berlekamp-Massey algorithm
# Adapted from www.cs.wisc.edu/~cs435-1/bermas.m
# Transliteration of
#   Massey, "Shift-Register Synthesis and BCH Decoding,"
#   IEEE Trans. Inform. Theory, 15(1):122-127, 1969.
# Input: P, either 0 or a prime
#           If P>0 then we work over the field K = Z/Z[P] (mod P)
#           else we work over the field K = Q (rationals)
#        N, a positive integer
#        s, a list of >= 2*N terms in K
#        x, a formal variable
# Returns: Unique monic annihilator of minimum degree, over K[x].

 BM := proc(s, N, P, x)
   local C,B,T,L,k,i,n,d,b,safemod;
   ASSERT(nops(s) = 2*N);
   safemod := (exp, P) -> `if`(P=0, exp, exp mod P);
   B := 1;
   C := 1;
   L := 0;
   k := 1;
   b := 1;
   for n from 0 to 2*N-1 do
     d := s[n+1];
     for i from 1 to L do
       d := safemod(d + coeff(C,x^i)*s[n-i+1], P);
     od;
     if d=0 then k := k+1 fi;
     if (d <> 0 and 2*L > n) then
       C := safemod(expand(C - d*x^k*B/b), P);
       k := k+1;
     fi;
     if (d <> 0 and 2*L <= n) then
       T := C;
       C := safemod(expand(C - d*x^k*B/b), P);
       B := T;
       L := n+1-L;
       k := 1;
       b := d;
     fi;
   od;
   return C;
 end:


# Compute lambda(z) usign the Euclidean algorithm
BMEA := proc(v::list,p::posint,z::name) local n,m,R0,R1,V0,V1,i;
  n := iquo( nops(v), 2 );
  m := 2*n-1;
  R0 := z^(2*n);
  R1 := add( v[m+1-i]*z^i, i=0..m );
lprint("R0=",R0);
lprint("R1=",R1);
  V0 := 0;
  V1 := 1;
  while n <= degree(R1,z) do
     R0,R1 := R1,Rem(R0,R1,z,'Q') mod p;
lprint("R0=",R0);
lprint("R1=",R1);
     V0,V1 := V1,Expand(V0-Q*V1) mod p;
lprint("V0=",V0);
lprint("V1=",V1);
  od;
  i := 1/lcoeff(V1,z) mod p;
  i*V1 mod p;
end:

testing := true:

if testing then

p := 103:
f := 3*x^3+5*x^6+7*x^11;
print("f=",f);
e := {3,6,11};
alpha := numtheory[primroot](p);# primitive root of Z_103
T := 4;
v := [seq( Eval(f,x=alpha^i) mod p, i=0..2*T-1 )];# list of numbers 
Lambda := BMEA(v,p,z);
R := Roots(Lambda) mod p;
[alpha^3, alpha^6, alpha^11 ] mod p;
L := { seq( numtheory[mlog]( r[1], alpha, p ), r in R ) };
if L = e then printf("Passed\n"); else printf("Failed\n"); fi;


#  P := 103:
#  d := 4:
#  num := 21+83*x+90*x^2+4*x^3: # degree < d
#  den := 1+11*x+23*x^2+58*x^3+69*x^4: # monic, degree <= d
#  f := series(num/den, x=0, 4*d) mod P;
#  s := [seq(coeff(f, x, i), i=0..13)];
#  g := BM(s, d, P, x);
#  h := BMEA(s,P,x);
#  Factor(h) mod P;

p := 997;
f := randpoly(x,degree=100,terms=5) mod p;
e := {seq( degree(t), t = f )};
f := unapply(f,x);
alpha := numtheory[primroot](p);
a := [seq( alpha^i mod p, i=0..19 )];
s := [seq( f(alpha^i mod p) mod p, i=0..19 )];
g := BM(s, nops(s)/2, p, x);
Lambda := BMEA(s,p,x);
R := Roots(Lambda) mod p;
R := [seq( r[1], r=R )];
L := {seq( numtheory[mlog](r,alpha,p), r=R )};
if L = e then printf("Passed\n"); else printf("Failed\n"); fi;

fi;

#F := unapply( add( c[i]*x^L[i], i=1..nops(R) ), x );
#n := nops(R);
#eqns := {seq( F(alpha^i mod p) mod p = s[i+1], i=0..n-1)};
#sols := msolve(eqns,p);
#eval( F(x), sols );
#f(x) mod p;

#
#roots=R;
#k := [seq( alpha^L[i] mod p, i=1..4 )]:
#'k'=k;
#M := mul(x-k[i],i=1..n):
#M := Expand(M) mod p:
#'M'=M;
#'L'=Lambda;
#
#P := Array(1..n):
#X := Vector(n):
#U := Vector(n);
#for j to n do
#    q := Quo(M,x-k[j],x) mod p;
#    print(degree(q));
#    U[j] := 1/eval(q,x=k[j]) mod p;
#    P[j] := q;
#    X[j] := U[j]*add( s[i]*coeff(P[j],x,i-1), i=1..n ) mod p;
#od:
#'X' = X;
#F(x);
#f(x);
#C := proc(i,j) coeff(P[i],x,j-1) end:
#M := Matrix(n,n,C);
#w := Vector(s[1..4]);
#c := M.w mod p:
#for i to n do c[i] := c[i]*U[i] mod p od:
#c=X;
#
#
##alpha := 7;
##p := 997;
##seq( alpha^i mod p, i=0..9 );
##f := (x-alpha^3)*(x-alpha^4) mod p;
##a := [seq( eval(f,x=alpha^i mod p) mod p, i=0..19 )];
##for i from 2 to 10 by 2 do Lambda := BMEA(s[1..i],p,x); Roots(Lambda) mod p; od;
#
#with(numtheory):
#f := x+x^10;
#p := 11;
#alpha := 2;
#a := [seq( eval(f,x=alpha^i mod p), i=0..5 )];
#Lambda := BMEA(a,p,x); 
#R := Roots(Lambda) mod p;
#seq( mlog(r[1],alpha,p), r=R );
#
