with( LinearAlgebra ):
LLL := proc( X, alpha )
local Y, Ystar, gstar, delta, mu, nu, xi,
i, j, k, l, m, n, r, t:
Y := copy( X ):
m := RowDimension( Y ):
n := ColumnDimension( Y ):
Ystar := Matrix( m, n ):
mu := Matrix( m, m ):
gstar := Vector( m ):
for i to m do
for t to n do Ystar[i,t] := Y[i,t] od:
for j to i-1 do
mu[i,j] :=
DotProduct( Row(Y,i), Row(Ystar,j) ) / gstar[j]:
Ystar := RowOperation( Ystar,[i,j], -mu[i,j] )
od:
gstar[i] := DotProduct( Row(Ystar,i), Row(Ystar,i) )
od:

k := 2:
while k <= m do
if abs(mu[k,k-1]) > 1/2 then
r := ceil( mu[k,k-1]-1/2 ):
Y := RowOperation( Y, [k,k-1], -r ):
for j to k-2 do mu[k,j] := mu[k,j] - r*mu[k-1,j] od:
mu[k,k-1] := mu[k,k-1] - r
fi:
if gstar[k] >= (alpha-mu[k,k-1]^2) * gstar[k-1] then
for l from k-2 to 1 by -1 do
if abs(mu[k,l]) > 1/2 then
r := ceil( mu[k,l]-1/2 ):
Y := RowOperation( Y, [k,l], -r ):
for j to l-1 do mu[k,j] := mu[k,j] - r*mu[l,j] od:
mu[k,l] := mu[k,l] - r
fi
od:
k := k + 1
else
nu := mu[k,k-1]:
delta := gstar[k] + nu^2*gstar[k-1]:
mu[k,k-1] := nu*gstar[k-1]/delta:
gstar[k] := gstar[k-1]*gstar[k]/delta:
gstar[k-1] := delta:
Y := RowOperation( Y, [k-1,k] ):
for j to k-2 do
t := mu[k-1,j]: mu[k-1,j] := mu[k,j]: mu[k,j] := t
od:
for i from k+1 to m do
xi := mu[i,k]:
mu[i,k] := mu[i,k-1] - nu*mu[i,k]:
mu[i,k-1] := mu[k,k-1]*mu[i,k] + xi
od:
if k > 2 then k := k - 1 fi
fi
od:
RETURN( Y )
end: