read "./TESTCASE_GENERATOR.mpl":
read "./0construct_BB.mpl":
read "./1_MRFI.mpl":
read "./2_NDSA.mpl":
read "./3_projection_phi.mpl":
read "./4_MQRFR.mpl":
read "./5_BMEA.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solve.mpl":
read "./8_construct_final_polynomial.mpl":



# test_case:="numerator_zero":
# THis is not working. Monomials are reconstructed properly but the coefficients are not. Not working for 25 variables. Try increasing prime. 
# num_var:=23:
# num_terms:=1432:
# den_terms:=1747:
# test_case:="rand":
# num_var:=10:
# num_terms:=17:
# den_terms:=11:
# vars,p,T,ff,gg:=data_generator(test_case,num_var,num_terms,den_terms):
test_case:=1:
vars,p,T,ff,gg:=data_generator(test_case):
num_var:=nops(vars):
print("vars=",vars):
print("p=",p):
print("T=",T):
print("ff=",ff mod p):
print("gg=",gg mod p):


B:=Construct_Rational_Blackbox(ff,gg,vars):
Numerator,Denominator:=MRFI(B,num_var,vars,p):
# print("ff= ",ff mod p):
# print("gg= ",gg):
# print("__________________________________________"):
r:=rand(p):
test_point:=[seq(r(),i=1..num_var)]:

print("checking "):
Eval(Numerator/Denominator ,{seq(vars[v]=test_point[v],v=1..num_var)}) mod p-B(test_point,p);

