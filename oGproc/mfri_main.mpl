read "./0_data_gen.mpl":
read "./1_construct_BB.mpl":
read "./get_primes.mpl":
read "./projection_phi.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
read "./4_BMEA.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solve.mpl":
read "./8_construct_final_polynomial.mpl":
read "./exp_vec.mpl":
read "./mfri_2.mpl":

test_case:="rand":
# test_case:="numerator_zero":
# test_case:=1:
num_var:=21:
num_terms:=160:
den_terms:=210:
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
