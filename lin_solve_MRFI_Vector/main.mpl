read "./0_construct_BB.mpl":
read "./1_MRFI.mpl":
read "./2a_NDSA.mpl":
read "./2b_Deterministic_NDSA.mpl":
read "./3_projection_phi.mpl":
read "./4_MQRFR.mpl":
read "./5_BMEA.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solver.mpl":
read "./8_construct_final_polynomial.mpl":
read "./9_rvr.mpl":
read "./data_gen.mpl":
read "./helpers.mpl":

p:= 2^31 - 1:

# Vars,F,G,num_vars,num_eqn,params := get_data(1):
# counter := 0:
# B:= Construct_Rational_Blackbox(F,G,Vars):



# test_case:="rand":
# num_var:=3:
# num_terms:=11:
# den_terms:=9:
# Vars,F,G,num_vars,num_eqn,params,field:=get_data(test_case,num_var,num_terms,den_terms):
# counter := 0:
# B:= Construct_Rational_Blackbox(F,G,Vars):


test_case:="rat_rand":
num_var:=3:
num_terms:=11:
den_terms:=9:
num_coeff_bound:=20:
den_coeff_bound:=20:
Vars,F,G,num_vars,num_eqn,params,field:=get_data(test_case,num_var,num_terms,den_terms,num_coeff_bound,den_coeff_bound):
counter := 0:
B:= Construct_Rational_Blackbox(F,G,Vars):

# lprint("Variables:", Vars):
# lprint("Numerator F:", F):
# lprint("Denominator G:", G):
# print("num_eqn =",num_eqn):



# test_case:="example":
# test_case:="small_sys_low_deg":
# test_case:="bspline":
# test_case:="small_Sys":

# Sys, Vars, params, num_vars, num_eqn,field:= get_data(test_case):
# counter := 0:
# B := Constuct_Sys_Blackbox(Sys, Vars, params):




# print("Number of equations:", num_eqn):
# print("Number of parameters:", num_vars):
# Create black box
try
NC,DC,NM,DM := MRFI(B, num_vars, num_eqn, params, p):
    catch:
    lprint("ERROR:", lasterror()):
end try:

NCQ:=table():
DCQ:=table():
for i from 1 to num_eqn do 
    NCQ[i]:=RVR(NC[i],p,field):
    DCQ[i]:=RVR(DC[i],p,field):
end do:
N_Fin:=table():
D_Fin:=table():
Rat_recon:=table():




print("======================================================"):
print("Displaying the results"):

if(num_eqn >1)then 
    og_soln:=get_eqn(Sys,Vars):
    og_unordered_soln:=convert(og_soln,list):
    og_soln:=reording(og_unordered_soln,nops(Sys)):
    fin_rat_recon:=Vector(convert(Rat_recon,list)):
    for i from 1 to num_eqn do 
        print("x",i,"="):
        N_Fin[i]:=construct_final_polynomial(NCQ[i],NM[i]):
        D_Fin[i]:=construct_final_polynomial(DCQ[i],DM[i]):
        Rat_recon[i]:=N_Fin[i]/D_Fin[i]:
        print("Rat_recon= ",Rat_recon[i]):
        print("original_soln =",og_soln[i]):
    end do:
    elif num_eqn =1 then 
        N_Fin[1]:=construct_final_polynomial(NCQ[1],NM[1]):
        D_Fin[1]:=construct_final_polynomial(DCQ[1],DM[1]):
        Rat_recon[1]:=N_Fin[1]/D_Fin[1]:
        MRFI_B:= Construct_Rational_Blackbox(N_Fin[1],D_Fin[1],Vars):
        print("Black box for interpolated rational polynomial: MFRI_B"):
        # fin_rat_recon:=Vector(convert(Rat_recon,list)):
        lprint("Rat_recon= ",sort(Rat_recon[1],Vars)):
    # print("lcoeff g = ",lcoeff(g,))
        lprint("Original equation =",sort(F/G,Vars)):
        R:=rand(p):
        testing_point:=[seq(R(),i=1..num_vars)]:
        print("testing_point = ",testing_point):
        print("B(testing_point,p)-MRFI_B(testing_point,p) = ", B(testing_point,p)-MRFI_B(testing_point,p)):
        print("field = ",field):
end if:

print("======================================================"):

lprint("Total Black Box Calls:", counter):