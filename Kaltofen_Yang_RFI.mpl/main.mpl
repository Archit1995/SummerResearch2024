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
read "./9_RCR.mpl":
read "./10_term_wise_images.mpl":


# test_case:="rand":
# test_case:="numerator_zero":
# test_case:="denominator_zero":
test_case:="rat_rand":
num_var:=7:
num_terms:=14:
den_terms:=11:

p1:=2^31-1:
p2:=2^31+11:
p3:=ithprime(9000000):
P:=[p1,p2,p3]:
# THis is not working. Monomials are reconstructed properly but the coefficients are not. Not working for 25 variables. Try increasing prime. 
# num_var:=23:
# num_terms:=1432:
# den_terms:=1747:
# num_var:=10:
# num_terms:=17:
# den_terms:=11:
vars,T,ff,gg:=data_generator(test_case,num_var,num_terms,den_terms,p1):
print("vars=",vars):
print("ff= ",ff):
print("gg= ",gg):

B:=Construct_Rational_Blackbox(ff,gg,vars):
print("B= ",B):
Numerator,Denominator:=RCR(B,num_var,vars,P):
print("Numerator= ",Numerator):
print("Denominator= ",Denominator):
print("Numerator-ff= ",Numerator-ff):
print("Denominator-gg= ",Denominator-gg):