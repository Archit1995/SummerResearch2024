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


test_case:="rand":
# test_case:="numerator_zero":
# test_case:=1:
num_var:=7:
num_terms:=10:
den_terms:=21:
vars,p,T,ff,gg:=data_generator(test_case,num_var,num_terms,den_terms):
print("vars=",vars):
print("p=",p):


B:=Construct_Rational_Blackbox(ff,gg,vars):
Numerator,Denominator:=MRFI(B,num_var,vars,p):
print("ff= ",ff mod p):
print("gg= ",gg):
print("__________________________________________"):
print("Numerator= ",Numerator-ff):
print("Denominator= ",Denominator-gg):
