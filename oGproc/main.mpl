read "./0_data_gen.mpl":
read "./1_construct_BB.mpl":
read "./projection_phi.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
read "./4_BMEA.mpl":
read "./5_get_lambda.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solve.mpl":
read "./8_construct_final_polynomial.mpl":
read "./mrfi.mpl":
# read "./get_sigma.mpl":
test_case:="denominator_zero":
# test_case:="numerator_zero":
# test_case:=1:
num_var,vars,p,T,ff,gg:=data_generator(test_case):
print("p=",p):
print("ff= ",ff mod p):
print("gg= ",gg):
# test:=8:
B:=Construct_Rational_Blackbox(ff,gg,vars):
# print("B: ",B([2,3,5],p)):
r_:=rand(p):
# shift_:=[r_()*r_() mod p]:
shift_:=[r_(),r_()]:
# shift_:=[r_()]:
# print("shift_: ",shift_):

anchor_point:=[2,3,5]:

# anchor_point:=[r_(),r_(),r_()]:
# Get degree of the polynomials and number of points
# We need sigma=[sig1,sig2,...sign], beta=[beta1,beta2,...beta_n-1]
# get_degree(B,sigma,beta,p) 
# terms_num,terms_den,lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:=MRFI(B,num_var,vars,p):
terms_num,terms_den,lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:=get_lambda(B,anchor_point,shift_,T,num_var,p,test_case):
Roots_num := [ seq(r[1], r in R_num ) ]:
Roots_den := [ seq(r[1], r in R_den ) ]:
 
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
num:= construct_final_polynomial(coeff_num,num_mono):
den:= construct_final_polynomial(coeff_den,den_mono):
# print("num: ",num):
# print("den: ",den):
den_lc:=lcoeff(den,order=grlex(x,y,z)):
den_lc_inv:=1/den_lc mod p:
print("den_lc: ",den_lc):
print("num= ",num*den_lc_inv mod p):
print("den= ",den*den_lc_inv mod p):



# sing primes [p1,....pn] where pn is the largest then make sure
#  p> pn^(total degree of polynomial we are interpolating)

# There could be trouble if if 2,3,5 is a root. 