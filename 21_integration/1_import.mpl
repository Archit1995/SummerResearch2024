read "/home/archit/summer2024_Research/13_final_sep/1_final_sep.mpl":
vars:={x,y}:
num_var:=numelems(vars):
p:=2^31-1:
# Test 1 
num_deg:=3:
denom_deg:=2:
ff:=randpoly(vars,sparse,degree=num_deg) mod p:
gg:=randpoly(vars,sparse,degree=denom_deg) mod p:
gg:=gg/lcoeff(gg) mod p:
print("ff ",ff):
print("gg ",gg):
# degree n=70 and d=52 CAUSING PROBLEM
B:=Construct_Blackbox(ff,gg,vars):
print(B):
f,g,lc_g,sigma_,shift_:=Early_termination_seperation(B,p):
Phi:=shift_[1]*x-shift_[1]*sigma_[1]+sigma_[2] mod p:
f_x:=Expand(subs(y=Phi,ff)) mod p:
print("f_x= ",f_x):
g_x:=Expand(subs(y=Phi,gg)) mod p:
print("f= ",f):
f_x:=f_x/lcoeff(g_x) mod p:
g_x:=g_x/lcoeff(g_x) mod p:
print("g_x= ",g_x):
# Phi:=shift_[1]*x-shift_[1]*sigma_[1]+sigma_[2] mod p:
# f_x:=Expand(subs(y=Phi,ff)) mod p:
# print("f_x= ",f_x):
# g_x:=Expand(subs(y=Phi,gg)) mod p:
# print("f= ",f):
# f_x:=f_x/lcoeff(g_x) mod p:
# g_x:=g_x/lcoeff(g_x) mod p:
print("checking numerator: ",f-f_x):
print("checking denominator: ",g-g_x):

