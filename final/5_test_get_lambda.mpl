read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
read "./4_BMEA.mpl":
read "./5_get_lambda.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solve.mpl":
p:=2^31-1:
# p:=19:
T:=4:
# p:=811:
# p:=19:
ff:=x*y+1:
gg:=x+y:
# ff:=x^3*y^3+x^2*y^2+x*y+1:
# gg:=x^2+y^2+x+y+3:
num_var:=2:
vars:={x,y}:
test:=8:
B:=Construct_Rational_Blackbox(ff,gg,vars):
anchor_point:=[2,3]:
shift_:=[10]:
terms_num,terms_den,lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:=get_lambda(B,anchor_point,shift_,T,num_var,p):
Roots_num := [ seq(r[1], r in R_num ) ]:
Roots_den := [ seq(r[1], r in R_den ) ]:
print("Roots_=",Roots_):
num_mono:= generate_monomials(Roots_num,num_var,anchor_point,vars):
print("num_mono: ",num_mono):
den_mono:= generate_monomials(Roots_den,num_var,anchor_point,vars):
print("den_mono: ",den_mono):
print("terms_num: ",terms_num):
# print(type(Roots_num,list)):
print("Roots_num: ",Roots_num):
coeff_num:= Zippel_Vandermonde_solver(num_vals,terms_num,Roots_num,lambda_num,p):
print("coeff_num: ",coeff_num):
coeff_den:= Zippel_Vandermonde_solver(den_vals,terms_den,Roots_den,lambda_den,p):
print("coeff_den: ",coeff_den):

# Ben Or Tiwari
# If u

# sing primes [p1,....pn] where pn is the largest then make sure
#  p> pn^(total degree of polynomial we are interpolating)

# There could be trouble if if 2,3,5 is a root. 