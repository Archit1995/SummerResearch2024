read "./0data_gen.mpl":
read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
with(NumberTheory):
vars:={x}:

p:=31:
T:=4:
ff:=x^2+5*x+1:
gg:=s*x+2;
actual_lc:=lcoeff(gg):
Eval(ff,x=2) mod p;
Eval(gg,x=2) mod p;
# ff,gg,actual_lc:=data_generator(0,deg,p):
# print("ff= ",ff):
# print("gg= ",gg):
# B:=Construct_Rational_Blackbox(ff,gg,vars):
# print("B= ",B):
# B(2,p);
# rat:=ff/gg:
# Eval(rat,x=2) mod p;
# f,g,mqrfr_lc:=evaluation_num_den(B,T,p):
# print("f= ",f):
# print("g= ",g):
# print("actual_lc= ",actual_lc):
# print("mqrfr_lc= ",mqrfr_lc):
# Question? check 1/lcoeff(num)