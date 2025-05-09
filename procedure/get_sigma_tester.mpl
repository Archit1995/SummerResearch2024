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
read "./get_sigma.mpl":
read "./generate_evaluation_primes.mpl":

test_case:="unlucky_primes":
# test_case:="numerator_zero":
# test_case:=1:
num_var,vars,p,T,ff,gg:=data_generator(test_case):
print("p=",p):
print("ff= ",ff):
print("gg= ",gg):
B:=Construct_Rational_Blackbox(ff,gg,vars):
print("trying to get sigma"):
anchor_point:=get_sigma(B,num_var,p):
print("anchor_point=",anchor_point):