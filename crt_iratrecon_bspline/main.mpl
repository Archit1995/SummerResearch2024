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

p1:=2^31-1:
p2:=nextprime(p1):
# p3:=nextprime(p2):


# Vars,F,G := get_data(1):
# B:= Construct_Rational_Blackbox(F,G,Vars):
# params := Vars:num_eqn:=nops(B([2,3], p)):
# num_vars := nops(params):
# lprint("Variables:", Vars):
# lprint("Numerator F:", F):
# lprint("Denominator G:", G):


# Sys, Vars, params, p := get_data("example"):
# Sys, Vars, params, p := get_data("small_Sys"):
# Sys, Vars, params, p := get_data("small_sys_low_deg"):
Sys, Vars, params := get_data("bspline"):
B := Constuct_Sys_Blackbox(Sys, Vars, params):
num_vars := nops(params):
num_eqn := nops(Vars):
lprint("Variables:", Vars):
lprint("Parameters:", params):
lprint("Prime:", p1):


# print("Number of equations:", num_eqn):
# print("Number of parameters:", num_vars):
# Create black box

p4:=nextprime(p3):
p5:=nextprime(p4):
Mp:=p1*p2:
# Mp:=p1*p2*p3*p4*p5:

try
N1,D1,num_mono,den_mono := MRFI_Vector(B, num_vars, num_eqn, params, p1):
N2,D2,num_mono,den_mono := MRFI_Vector(B, num_vars, num_eqn, params, p2):
# N3,D3,num_mono,den_mono := MRFI_Vector(B, num_vars, num_eqn, params, p3):
# N4,D4,num_mono,den_mono := MRFI_Vector(B, num_vars, num_eqn, params, p4):
# N5,D5,num_mono,den_mono := MRFI_Vector(B, num_vars, num_eqn, params, p5):
print("N1= ",N1 );
# print("D1= ",D1 ):
# print("N2= ",N2 ):
# print("D2= ",D2 ):
N_crt:=table():
D_crt:=table():
for k from 1 to num_eqn do
    N_crt[k] := chrem([N1[k], N2[k]], [p1, p2]):
    D_crt[k] := chrem([D1[k], D2[k]], [p1, p2]):
end do:

# for k from 1 to num_eqn do
#     N_crt[k] := chrem([N1[k], N2[k], N3[k], N4[k]], [p1, p2,p3,p4]):
#     D_crt[k] := chrem([D1[k], D2[k], D3[k], D4[k]], [p1, p2,p3,p4]):
# end do:



# for k from 1 to num_eqn do
#     N_crt[k] := chrem([N1[k], N2[k], N3[k], N4[k], N5[k]], [p1, p2,p3,p4,p5]):
#     D_crt[k] := chrem([D1[k], D2[k], D3[k], D4[k], D5[k]], [p1, p2,p3,p4,p5]):
# end do:

# print("N_crt= ",N_crt ):
# print("D_crt= ",D_crt ):
# N_rat:=table():
# D_rat:=table():
# for k from 1 to num_eqn do
#     N_rat[k]:=iratrecon(N_crt[k],Mp):
#     if(N_rat[k]=FAIL)then N_rat_done[k]:=false; else N_rat_done[k]:=true; end if:
#     D_rat[k]:=iratrecon(D_crt[k],Mp):
#     if(D_rat[k]=FAIL)then D_rat_done[k]:=false; else D_rat_done[k]:=true; end if:
# end do:
# print("N_rat= ",N_rat ):
# print("D_rat= ",D_rat ):
U_crt:=table():
for k from 1 to num_eqn do
    inv_d := (D_crt[k]^(-1)) mod Mp;
    U_crt[k] := (N_crt[k] *inv_d) mod Mp;
end do;
print("U_crt= ",U_crt ):
result:=table():
for k from 1 to num_eqn do
    result[k]:=iratrecon(U_crt[k],Mp):
end do;
print("Result from iratrecon on U_crt= ",result ):
# Numerators:=table():
# Denominators:=table():
# for k from 1 to num_eqn do
#     if(N_rat[k]<>FAIL)then 
#         Numerators[k]:= construct_final_polynomial(N_rat[k], num_mono[k]):
#         print("Numerator", k, "before normalization:", Numerators[k] ):
#     end if:
#     if(D_rat[k]<>FAIL)then
#         Denominators[k]:= construct_final_polynomial(D_rat[k], den_mono[k]):
#     end if;
# end do:
# print("Numerators before normalization:,",Numerators ):
# print("Denominators before normalization:,",Denominators ):
# den_lc_inv:=table():
# for k from 1 to num_eqn do
#     if(D_rat[k]<>FAIL)then 
#         den_lc:=lcoeff(Denominators[k],order=grlex(seq(params[i],i=1..num_vars))):
#         print("den_lc for equation", k, "=", den_lc ):
#         den_lc_inv[k]:=1/den_lc;
#     end if;
# end do:
# print("den_lc_inv= ",den_lc_inv ):
# for k from 1 to num_eqn do
#     if(N_rat_done[k]=true and D_rat_done[k]=true)then
#         print("Normalizing equation", k ):
#         Numerators[k]:=Numerators[k]*den_lc_inv[k];
#         Denominators[k]:=Denominators[k]*den_lc_inv[k];
#     end if;
# end do:
# print("Recovered Rational Function Components:"):
# for k from 1 to num_eqn do
#     lprint("Component", k, ":"):
#     lprint("  Numerator:", Numerators[k]):
#     lprint("  Denominator:", Denominators[k]):
# end do:

    catch:
    lprint("ERROR:", lasterror()):
end try:
