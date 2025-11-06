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
read "./recover.mpl":



# test_case:="numerator_zero":
# THis is not working. Monomials are reconstructed properly but the coefficients are not. Not working for 25 variables. Try increasing prime. 
# num_var:=23:
# num_terms:=1432:
# den_terms:=1747:
# test_case:="rand":
# num_var:=10:


# num_terms:=17:
# den_terms:=11:

# num_terms:=11:
# den_terms:=9:

# num_var:=3:
# test_case:="rand":
# vars,ff,gg:=data_generator(test_case,num_var,num_terms,den_terms):
# print("vars=",vars):
# print("ff=",ff):
# print("gg=",gg):


test_case:=3:


vars,ff,gg:=data_generator(test_case):
num_var:=nops(vars):
print("vars=",vars):
print("ff=",ff):
print("gg=",gg):

T:=4:
p:=2^31-1:
counter:=0:
B:=Construct_Rational_Blackbox(ff,gg,vars):
N1,D1,num_mono,den_mono:=MRFI(B,num_var,vars,p):
print("N1= ",N1 ):
print("D1= ",D1 ):
num:=construct_final_polynomial(N1,num_mono):
print("num = ",num ):
den:=construct_final_polynomial(D1,den_mono):
print("den = ",den ):

print("_______________________________________________________________________________________"):


Nrvr:=RVR(N1,p):
print("Nrvr= ",Nrvr):
Drvr:=RVR(D1,p):
print("Drvr= ",Drvr):

num:=construct_final_polynomial(Nrvr,num_mono):
print("num = ",num ):
print("ff= ",ff):
den:=construct_final_polynomial(Drvr,den_mono):
print("den = ",den ):

print("gg= ",gg):
print("========================================================================================"):
og_rat:=ff/gg:
print("Original rational function = ",og_rat):
recovered_rat:=num/den:
print("Recovered rational function = ",recovered_rat):
if simplify(og_rat - recovered_rat) = 0 then
    print("Success: The original and recovered rational functions are identical."):
else
    print("Failure: The original and recovered rational functions differ."):
end if:
print("Counter value (number of blackbox calls): ", counter):