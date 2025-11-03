read "./0_construct_BB.mpl":
read "./1_MRFI.mpl":
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



# Sys, Vars, params, p := get_data("example"):
Sys, Vars, params, p := get_data("small_Sys"):
# Sys, Vars, params, p := get_data("small_sys_low_deg"):
# Sys, Vars, params, p := get_data("bspline"):
# Sys, Vars, params, p := get_data("vector_example"):
print("System"):
lprint("Eq 1=:", Sys[1]):
lprint("Eq 2=:", Sys[2]):
lprint("Variables:", Vars):
lprint("Parameters:", params):
lprint("Prime:", p):
num_vars := nops(params):
num_eqn := nops(Vars):
print("Number of equations:", num_eqn):
print("Number of parameters:", num_vars):
# Create black box
B := Constuct_Sys_Blackbox(Sys, Vars, params):
# B:= Constuct_Vector_Blackbox(Sys, Vars, params):
# print("B[2,3] =",B([2,3], p)):


try
Numerators, Denominators := MRFI_Vector(B, num_vars, num_eqn, params, p):

    for k from 1 to num_eqn do
        lprint(""):
        lprint("Component", k, ":"):
        lprint("  Numerator:", Numerators[k]):
        lprint("  Denominator:", Denominators[k]):
    end do:

    catch:
    lprint("ERROR:", lasterror()):
end try:
