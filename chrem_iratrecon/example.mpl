
p1:=2^31-1:
# p2:=2^31+11:
p2:=ithprime(1000):
P:=[p1,p2]:
# a:=3/4:
a:=7/19:
b:=[]:
r:=rand(p1):
c:=r():
b:=[op(b),seq(a mod P[i],i=1..nops(P))]:
b2:=1610612736 mod p2:
b:=[1610612736,b2]:
print("b=",b):
m:=product(P[i],i=1..nops(P)):
print("m=",m):
u:=chrem(b,P):
iratrecon(u,m);