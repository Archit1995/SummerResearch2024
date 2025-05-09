read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":

# f:=
# g:=x+2:
# B:=f/g:
p:=19:
# p:=2^31-1:
alpha:=modp([0,1,2,3,4],p):
print("alpha= ",alpha): 
# y:=[seq(eval(B,x=alpha[i])mod p,i=1..5)]:
y:=[8,1,9,4,2];
print("y= ",y):
m:=expand(x*(x-1)*(x-2)*(x-3)*(x-4)) mod p:
print("m= ",m):
u:=Interp(alpha,y,x)mod p:
# u:=u/lcoeff(u) mod p:
print("u= ",u):
num,den,q,g_lcoeff:=MQRFR(m,u,0,1,p):
print("num= ",num):
print("den= ",den):
print("q= ",q):
print("g_lcoeff= ",g_lcoeff):
gcd(num,den):
print("gcd(num,den)= ",gcd(num,den)):
