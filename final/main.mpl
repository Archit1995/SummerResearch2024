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
# num_var:=3:
# num_terms:=17:
# den_terms:=11:
# num_terms:=3:
# den_terms:=2:
# vars,ff,gg:=data_generator(test_case,num_var,num_terms,den_terms):
# test_case:="bspline_small_sys_low_deg2":
# test_case:="bspline_small_sys_low_deg3":
test_case:=3:
# test_case:="rand":
# p1:=2^31-1:
# p2:=nextprime(2^31-1):
# Mp:=p1*p2:

T:=4:
vars,ff,gg:=data_generator(test_case):
num_var:=nops(vars):
print("vars=",vars):
print("ff=",ff):
print("gg=",gg):


p1:=2^31-1:

# print("B([2,3], p1) =",B([2,3,5], p1)):

# # primes_to_use := [2147483647, 2147483659, 2147483693, 2147483699, 2147483701]:
# num_primes := 5:
# primes_to_use:=[seq(i, i=0..num_primes-1)]:
# primes_to_use[1]:=p1:
# for i from 2 to num_primes do
#     pnext:=nextprime(primes_to_use[i-1]):
#     primes_to_use[i]:=pnext:
# end do:
# print("Primes to use =", primes_to_use):
# num_primes := nops(primes_to_use):
# N_list := []:
# D_list := []:
# prime_list := []:
# for i from 1 to num_primes do
#     p := primes_to_use[i]:
#     print("========================================================"):
#     print("i=", i):
#     print("Using prime p =", p):
#     B:=Construct_Rational_Blackbox(ff,gg,vars):
#     N_i, D_i, num_mono, den_mono := MRFI(B, num_var, vars, p):
#     print("N_i = ", N_i):
#     print("D_i = ", D_i):
#     N_list := [op(N_list), N_i]:
#     D_list := [op(D_list), D_i]:
#     prime_list := [op(prime_list), p]:
# end do:
# print("N_list = ", N_list):
# print("D_list = ", D_list):
# print("prime_list = ", prime_list):

# print(""):
# print("========================================"):
# print("ALL RUNS COMPLETE"):
# print("========================================"):
# print("Now performing rational reconstruction...");
# `mod`  := mods;
# comnined_N:=chrem([op(N_list)], prime_list):
# comnined_D:=chrem([op(D_list)], prime_list):
# print("Combined Numerator =", comnined_N):
# print("Combined Denominator =", comnined_D):

# # Calculate combined modulus
# M := mul(p, p in primes_to_use):
# print("Combined modulus M =", M):
# print("Decimal digits:", length(convert(M, string))):

# print(""):
# print("========================================"):
# print("RATIONAL RECONSTRUCTION"):
# print("========================================"):
# # final_numerator := iratrecon(comnined_N, M):
# final_numerator :=[seq(iratrecon(comnined_N[i], M), i=1..nops(comnined_N))]:
# # final_denominator := iratrecon(comnined_D, M):
# final_denominator :=[seq(iratrecon(comnined_D[i], M), i=1..nops(comnined_D))]:
# print("Final Numerator Coefficients =", final_numerator):
# print("Final Denominator Coefficients =", final_denominator):



# num := construct_final_polynomial(final_numerator, num_mono):
# den := construct_final_polynomial(final_denominator, den_mono):
# den_lc := lcoeff(den, order = grlex(seq(vars[i], i = 1 .. num_var))):
# inv_den_lc := 1 / den_lc:
# num := num * inv_den_lc ;
# den := den * inv_den_lc ;
# print("Final Numerator = ", num):
# print("Final Denominator = ", den):
# gg_lc := lcoeff(gg, order = grlex(seq(vars[i], i = 1 .. num_var))):
# inv_gg_lc := 1 / gg_lc:
# print("ff= ", ff*inv_gg_lc):
# print("gg= ", gg*inv_gg_lc):
# print("__________________________________________"):
B:=Construct_Rational_Blackbox(ff,gg,vars):
N1,D1,num_mono,den_mono:=MRFI(B,num_var,vars,p1):
print("N1= ",N1 ):
print("D1= ",D1 ):
print("_______________________________________________________________________________________"):
p2:=nextprime(p1):
N2,D2,num_mono,den_mono:=MRFI(B,num_var,vars,p2):
print("N2= ",N2 ):
print("D2= ",D2 ):
print("[N1 N2]= ",[N1,N2] ):
print("[D1 D2]= ",[D1,D2] ):
print("_______________________________________________________________________________________"):
`mod`  := mods;

N_crt := chrem([N1, N2], [p1, p2]):
print("N_crt= ",N_crt ):
D_crt := chrem([D1, D2], [p1, p2]):
print("D_crt= ",D_crt ):
print("_______________________________________________________________________________________"):
Mp:=mul(p, p in [p1,p2]):
print("Mp= ",Mp ):
Nrat:=iratrecon(N_crt,Mp);
Drat:=iratrecon(D_crt,Mp);
print("Nrat= ",Nrat ):
print("Drat= ",Drat ):
num:= construct_final_polynomial(Nrat,num_mono):
den:= construct_final_polynomial(Drat,den_mono):


den_lc:=lcoeff(den,order=grlex(seq(vars[i],i=1..num_var))):
inv_den_lc:=1/den_lc;


num:=num*inv_den_lc ;
den:=den*inv_den_lc ;
print("num = ",num ):
print("den = ",den ):


print("ff= ",ff):
print("gg= ",gg):
# print("__________________________________________"):
# r:=rand(p):
# test_point:=[seq(r(),i=1..num_var)]:

# print("checking "):
# Eval(Numerator/Denominator ,{seq(vars[v]=test_point[v],v=1..num_var)}) mod p-B(test_point,p);

