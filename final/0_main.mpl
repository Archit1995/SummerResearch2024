read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":
read "./2_seperation.mpl":
num_var:=2:
vars:={x,y}:
p:=2^31-1:
num_deg:=4:
denom_deg:=3:
ff:=randpoly(vars,sparse,degree=num_deg) mod p:
gg:=randpoly(vars,sparse,degree=denom_deg) mod p:
print("ff ",ff):
print("gg ",gg):
B:=Construct_Rational_Blackbox(ff,gg,vars);
test_separation(ff,gg,B,p):
# f,g,lc_g,sigma_,shift_:=Early_termination_seperation(B,p):
# f_x:=f_x/lcoeff(g_x) mod p:
# print("f_x= ",f_x):
# g_x:=g_x/lcoeff(g_x) mod p:
# # print("lc_g= ",lc_g mod p):
# print("g_x= ",g_x):
# # print("g= ",g):
# print("checking numerator: ",f-f_x):
# print("checking denominator: ",g-g_x):

