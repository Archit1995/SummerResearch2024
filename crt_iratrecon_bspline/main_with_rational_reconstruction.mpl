########################################################
# MAIN WITH MULTI-PRIME RATIONAL RECONSTRUCTION
########################################################

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

########################################################
# RATIONAL RECONSTRUCTION PROCEDURE
########################################################

Reconstruct_Rational_Poly := proc(poly_images::list, monomial_list::list, 
                                   primes::list, vars::list)
    local M, num_monoms, i, coeffs_at_prime, coeff_images, j,
          combined, rat_coeff, reconstructed_coeffs, final_poly:
    
    M := mul(p, p in primes):
    num_monoms := nops(monomial_list):
    
    # Extract coefficients from each polynomial image
    coeff_images := []:
    for i from 1 to nops(poly_images) do
        coeffs_at_prime := []:
        for j from 1 to num_monoms do
            # Get coefficient of j-th monomial in i-th polynomial
            coeffs_at_prime := [op(coeffs_at_prime), 
                               coeff(poly_images[i], monomial_list[j])]:
        end do:
        coeff_images := [op(coeff_images), coeffs_at_prime]:
    end do:
    
    # Reconstruct each coefficient
    reconstructed_coeffs := []:
    for j from 1 to num_monoms do
        # Get j-th coefficient from all primes
        combined := chrem([seq(coeff_images[i][j], i=1..nops(primes))], primes):
        
        # Try rational reconstruction
        rat_coeff := iratrecon(combined, M):
        
        if rat_coeff = FAIL then
            print("WARNING: Coefficient", j, "reconstruction failed"):
            print("  Using ratrecon with bounds..."):
            rat_coeff := ratrecon(combined, M, 10^10, 10^10):
            if rat_coeff = FAIL then
                print("  Still failed. Using integer value."):
                rat_coeff := combined:
            end if:
        end if:
        
        reconstructed_coeffs := [op(reconstructed_coeffs), rat_coeff]:
    end do:
    
    # Build final polynomial
    final_poly := add(reconstructed_coeffs[i] * monomial_list[i], 
                     i = 1..num_monoms):
    
    return final_poly:
end proc:

########################################################
# MAIN EXECUTION WITH MULTIPLE PRIMES
########################################################

# Choose test case
Sys, Vars, params, p_base := get_data("vector_example"):

print("System:"):
lprint("Eq 1:", Sys[1]):
lprint("Eq 2:", Sys[2]):
lprint("Variables:", Vars):
lprint("Parameters:", params):

num_vars := nops(params):
num_eqn := nops(Vars):
print("Number of equations:", num_eqn):
print("Number of parameters:", num_vars):

# Define multiple primes
primes_to_use := [2147483647, 2147483659, 2147483693, 2147483699, 2147483701]:
num_primes := nops(primes_to_use):

print(""):
print("========================================"):
print("MULTI-PRIME EXECUTION"):
print("Using", num_primes, "primes"):
print("Primes:", primes_to_use):
print("========================================"):

# Storage for results from each prime
all_numerators := []:
all_denominators := []:

# Run MRFI for each prime
for i from 1 to num_primes do
    p := primes_to_use[i]:
    
    print(""):
    print("========================================"):
    print("RUN", i, "OF", num_primes):
    print("Prime p =", p):
    print("========================================"):
    
    # Create black box for this prime
    B := Constuct_Vector_Blackbox(Sys, Vars, params):
    
    try
        Numerators, Denominators := MRFI_Vector(B, num_vars, num_eqn, params, p):
        
        print(""):
        print("Run", i, "completed successfully"):
        
        # Store results
        all_numerators := [op(all_numerators), Numerators]:
        all_denominators := [op(all_denominators), Denominators]:
        
    catch:
        error "MRFI failed for prime", p, ":", lasterror():
    end try:
end do:

print(""):
print("========================================"):
print("ALL RUNS COMPLETE"):
print("========================================"):
print("Now performing rational reconstruction...");

# Calculate combined modulus
M := mul(p, p in primes_to_use):
print("Combined modulus M =", M):
print("Decimal digits:", length(convert(M, string))):

# Reconstruct each component
print(""):
print("========================================"):
print("RATIONAL RECONSTRUCTION"):
print("========================================"):

final_numerators := []:
final_denominators := []:

for k from 1 to num_eqn do
    print(""):
    print("Component", k, ":"):
    
    # Get numerator polynomials for this component across all primes
    num_polys := [seq(all_numerators[i][k], i=1..num_primes)]:
    den_polys := [seq(all_denominators[i][k], i=1..num_primes)]:
    
    # Extract monomials from first polynomial (should be same across all)
    num_terms := [coeffs(num_polys[1], params, 'num_monoms')]:
    num_monoms := [num_monoms]:
    
    den_terms := [coeffs(den_polys[1], params, 'den_monoms')]:
    den_monoms := [den_monoms]:
    
    print("  Numerator has", nops(num_monoms), "terms"):
    print("  Denominator has", nops(den_monoms), "terms"):
    
    # Reconstruct numerator
    final_num := Reconstruct_Rational_Poly(num_polys, num_monoms, 
                                           primes_to_use, params):
    
    # Reconstruct denominator  
    final_den := Reconstruct_Rational_Poly(den_polys, den_monoms,
                                           primes_to_use, params):
    
    final_numerators := [op(final_numerators), final_num]:
    final_denominators := [op(final_denominators), final_den]:
    
    print("  Reconstructed numerator:", final_num):
    print("  Reconstructed denominator:", final_den):
end do:

print(""):
print("========================================"):
print("FINAL RESULTS (with rational coefficients):"):
print("========================================"):

for k from 1 to num_eqn do
    print(""):
    print("Component", k, ":"):
    print("  Numerator:", final_numerators[k]):
    print("  Denominator:", final_denominators[k]):
    print("  Rational function:", final_numerators[k] / final_denominators[k]):
end do:
