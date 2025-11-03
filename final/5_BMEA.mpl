BMEA := proc(v::list,p::posint,Z::name) # Might need correction because roots of lambda polynomial are not factors of 2,3,5 for all degrees
  local n,m,r,t,q,j:
  print("In BMEA");
  print("v=",v);    
  n := iquo( nops(v), 2 ):
  print("n=",n);
  m := 2*n-1:
  print("m=",m);
  r[0] := Z^(2*n):
  print("r0=",r[0]);
  r[1] := add( v[m+1-i]*Z^i, i=0..m ) mod p:
  print("r1=",r[1]);
# lprint("R0=",R0):
# lprint("R1=",R1):
  t[0] := 0:
  t[1] := 1:
  lprint("t0=",t[0]):
lprint("t1=",t[1]):
i:=1;
  while n <= degree(r[i],Z) do
  print("--------------------------------------"):
  print("Degree of r[",i,"]=",degree(r[i],Z)):
  # q[i]:= Quo(r[i-1],r[i],Z,'r[i+1]') mod p:
  r[i+1] := Rem(r[i-1],r[i],Z,'q[i]') mod p:
  lprint("r[",i-1,"]=",r[i-1]):
  lprint("r[",i,"]=",r[i]):
  lprint("r[",i+1,"]=",r[i+1]):
  lprint("q[",i,"]=",q[i]):
  # t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
  t[i+1] := Expand(t[i-1]-q[i]*t[i]) mod p:
  lprint("t[",i,"]=",t[i]):
  lprint("t[",i+1,"]=",t[i+1]):
  i:=i+1:
  end do:
  # lprint("q=",q):
  j := 1/lcoeff(t[i],Z) mod p:
  print("j=",j):
  return j*t[i] mod p:
end:

# BMEA := proc(v::list,p::posint,Z::name) # Might need correction because roots of lambda polynomial are not factors of 2,3,5 for all degrees
#   local n,m,R0,R1,V0,V1,i:
#   print("In BMEA");
#   print("v=",v);    
#   n := iquo( nops(v), 2 ):
#   print("n=",n);
#   m := 2*n-1:
#   print("m=",m);
#   R0 := Z^(2*n):
#   print("r0=",R0);
#   R1 := add( v[m+1-i]*Z^i, i=0..m ) mod p:
#   print("r1=",R1);
# # lprint("R0=",R0):
# # lprint("R1=",R1):
#   V0 := 0:
#   V1 := 1:
#   lprint("t0=",V0):
# lprint("t1=",V1):
#   while n <= degree(R1,Z) do
#      R0,R1 := R1,Rem(R0,R1,Z,'Q') mod p:
#     lprint("q=",Q):
# lprint("r0=",R0):
# lprint("r1=",R1):
#      V0,V1 := V1,Expand(V0-Q*V1) mod p:
# lprint("t0=",V0):
# lprint("t1=",V1):
#   od:
#   lprint("q=",Q):
#   i := 1/lcoeff(V1,Z) mod p:
#   print("i=",i):
#   return i*V1 mod p:
# end: