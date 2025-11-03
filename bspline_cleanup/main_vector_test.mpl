# Test file for vector MRFI using your modified files
read "./0_constructBB.mpl":
read "./1_MRFI_vector_complete.mpl":  # Using the cleaned up version
read "./2_NDSA.mpl":
read "./3_projection_phi.mpl":
read "./4_MQRFR.mpl":
read "./5_BMEA.mpl":
read "./6_generate_monomials.mpl":
read "./7_zippel_vandermonde_solve.mpl":
read "./8_construct_final_polynomial.mpl":
read "./reordering.mpl":
read "./data_gen.mpl":
read "./get_u.mpl":

print("==============================================="):
print("Testing Vector MRFI with Linear System Solution"):
print("==============================================="):

# Test Case 1: Your example system
print(""):
print("TEST CASE 1: Simple 2x2 system"):
print("---------------------------------"):

# Define the system that produces your example rational functions
# [x1 = -(y1^2-2*y1+3)/(y1*y2-1), x2 = -(y1*y2-y1-3*y2+1)/(y1*y2-1)]
Sys := {(y1*y2-1)*x1 + (y1^2-2*y1+3), (y1*y2-1)*x2 + (y1*y2-y1-3*y2+1)}:
Vars := {x1, x2}:
params := [y1, y2]:
p := 107:

print("System:", Sys):
print("Variables:", Vars):
print("Parameters:", params):

# Create the linear system black box
LBB := Constuct_Sys_Blackbox(Sys, Vars, params):

# Test the black box
test_point := [2, 3]:
result := LBB(test_point, p):
print("Black box evaluation at", test_point, "=", result):

# Run MRFI
num_vars := nops(params):
num_eqn := nops(Vars):

print(""):
print("Running MRFI_Vector..."):
print(""):

try
    Numerators, Denominators := MRFI_Vector(LBB, num_vars, num_eqn, params, p):
    
    print(""):
    print("RECOVERY COMPLETE!"):
    print("=================="):
    
    for k from 1 to num_eqn do
        print("Component", k, ":");
        print("  Recovered Numerator:", Numerators[k]):
        print("  Recovered Denominator:", Denominators[k]):
        print(""):
    end do:
    
    # Verify the recovery
    print("VERIFICATION:"):
    print("-------------"):
    
    # Expected results (from your example)
    expected_num1 := -(y1^2-2*y1+3) mod p:
    expected_num2 := -(y1*y2-y1-3*y2+1) mod p:
    expected_den := y1*y2-1 mod p:
    
    print("Expected:");
    print("  Numerator 1:", expected_num1):
    print("  Numerator 2:", expected_num2):
    print("  Denominator:", expected_den):
    
    # Test with random points
    r := rand(p):
    for test_num from 1 to 3 do
        test_pt := [r(), r()]:
        
        # Skip if denominator is zero
        if Eval(expected_den, {y1=test_pt[1], y2=test_pt[2]}) mod p = 0 then
            print("Skipping test point", test_pt, "(denominator zero)"):
            next:
        end if:
        
        print(""):
        print("Test point", test_num, ":", test_pt):
        
        # Evaluate original black box
        original_vals := LBB(test_pt, p):
        print("  Black box values:", original_vals):
        
        # Evaluate recovered functions
        recovered_vals := []:
        for k from 1 to num_eqn do
            num_eval := Eval(Numerators[k], {y1=test_pt[1], y2=test_pt[2]}) mod p:
            den_eval := Eval(Denominators[k], {y1=test_pt[1], y2=test_pt[2]}) mod p:
            if den_eval = 0 then
                print("  Recovered denominator zero at test point!"):
                break:
            end if:
            recovered_vals := [op(recovered_vals), num_eval/den_eval mod p]:
        end do:
        print("  Recovered values:", recovered_vals):
        
        # Check difference
        if nops(recovered_vals) = num_eqn then
            diff := [seq(original_vals[k] - recovered_vals[k] mod p, k=1..num_eqn)]:
            print("  Difference:", diff):
            
            if convert(diff, set) = {0} then
                print("  ✓ MATCH"):
            else
                print("  ✗ MISMATCH"):
            end if:
        end if:
    end do:
    
catch:
    print("ERROR in MRFI_Vector:", lasterror()):
end try:

# Test Case 2: Using your predefined test case
print(""):
print("==============================================="):
print("TEST CASE 2: Small system from data_gen"):
print("==============================================="):

Sys2, Vars2, params2, p2 := get_data("small_Sys"):
print("System:", Sys2):
print("Variables:", Vars2):
print("Parameters:", params2):

LBB2 := Constuct_Sys_Blackbox(Sys2, Vars2, params2):
num_vars2 := nops(params2):
num_eqn2 := nops(Vars2):

print(""):
print("Running MRFI_Vector on small_Sys..."):
print(""):

try
    Numerators2, Denominators2 := MRFI_Vector(LBB2, num_vars2, num_eqn2, params2, p2):
    
    print(""):
    print("RECOVERY COMPLETE for small_Sys!"):
    print("================================="):
    
    for k from 1 to num_eqn2 do
        print("Component", k, "(", convert(Vars2,list)[k], "):");
        print("  Numerator:", Numerators2[k]):
        print("  Denominator:", Denominators2[k]):
        print(""):
    end do:
    
catch:
    print("ERROR in MRFI_Vector for small_Sys:", lasterror()):
end try:

print(""):
print("==============================================="):
print("All tests complete"):
print("==============================================="):
