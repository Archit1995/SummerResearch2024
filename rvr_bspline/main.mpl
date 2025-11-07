read "./0_construct_BB.mpl":
# read "./1_MRFI.mpl":
read "./1_MRFI_hashtable.mpl":
read "./2_NDSA.mpl":
read "./Deterministic_NDSA.mpl":
read "./3_projection_phi.mpl":
read "./4_MQRFR.mpl":
read "./5_BMEA.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solver.mpl":
read "./8_construct_final_polynomial.mpl":
read "./data_gen.mpl":
read "./helpers.mpl":
read "./rvr.mpl":

p:= 2^31 - 1:

# Vars,F,G := get_data(1):
# counter := 0:
# B:= Construct_Rational_Blackbox(F,G,Vars):
# params := Vars:num_eqn:=nops(B([12,13,15], p)):
# num_vars := nops(params):
# lprint("Variables:", Vars):
# lprint("Numerator F:", F):
# lprint("Denominator G:", G):


# Sys, Vars, params := get_data("example"):

# Sys, Vars, params := get_data("small_sys_low_deg"):


# Sys, Vars, params := get_data("bspline"):

Sys, Vars, params:= get_data("small_Sys"):
counter := 0:
B := Constuct_Sys_Blackbox(Sys, Vars, params):
num_vars := nops(params):
num_eqn := nops(Vars):
lprint("Variables:", Vars):
lprint("Parameters:", params):



# print("Number of equations:", num_eqn):
# print("Number of parameters:", num_vars):
# Create black box




try
NC,DC,NM,DM := MRFI(B, num_vars, num_eqn, params, p):

    # for k from 1 to num_eqn do
    #     lprint(""):
    #     lprint("Component", k, ":"):
    #     lprint("  Numerator:", Numerators[k]):
    #     lprint("  Denominator:", Denominators[k]):
    # end do:

    catch:
    lprint("ERROR:", lasterror()):
end try:
NCQ:=table():
DCQ:=table():
for i from 1 to num_eqn do 
    NCQ[i]:=RVR(NC[i],p):
    DCQ[i]:=RVR(DC[i],p):
end do:
N_Fin:=table():
D_Fin:=table():
Rat_recon:=table():
for i from 1 to num_eqn do 
    N_Fin[i]:=construct_final_polynomial(NCQ[i],NM[i]):
    D_Fin[i]:=construct_final_polynomial(DCQ[i],DM[i]):
    Rat_recon[i]:=N_Fin[i]/D_Fin[i]:
end do:
# print(N_Fin):
# print(D_Fin):
fin_rat_recon:=Vector(convert(Rat_recon,list)):
for i from 1 to num_eqn do 
    print("x",i,"="):
    lprint(fin_rat_recon[i]):
end do:

lprint("Total Black Box Calls:", counter):