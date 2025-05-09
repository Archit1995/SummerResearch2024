BMEA := proc(v::list,p::posint,Z::name) # Might need correction because roots of lambda polynomial are not factors of 2,3,5 for all degrees
  local n,m,R0,R1,V0,V1,i:
  print("In BMEA");
  print("v=",v);    
  n := iquo( nops(v), 2 ):
  print("n=",n);
  m := 2*n-1:
  print("m=",m);
  R0 := Z^(2*n):
  print("R0=",R0);
  R1 := add( v[m+1-i]*Z^i, i=0..m ):
  print("R1=",R1);
# lprint("R0=",R0):
# lprint("R1=",R1):
  V0 := 0:
  V1 := 1:
  while n <= degree(R1,Z) do
     R0,R1 := R1,Rem(R0,R1,Z,'Q') mod p:
# lprint("R0=",R0):
# lprint("R1=",R1):
     V0,V1 := V1,Expand(V0-Q*V1) mod p:
# lprint("V0=",V0):
# lprint("V1=",V1):
  od:
  i := 1/lcoeff(V1,Z) mod p:
  return i*V1 mod p:
end: