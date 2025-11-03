# Test script for Complete MRFI Vector Package
read "./COMPLETE_MRFI_VECTOR_PACKAGE.mpl":

print(""):
print("========================================"):
print("TESTING MRFI FOR VECTOR BLACK BOXES"):
print("========================================"):
print(""):

# Test Case: Your example
# [x1 = -(y1^2-2*y1+3)/(y1*y2-1), x2 = -(y1*y2-y1-3*y2+1)/(y1*y2-1)]
print("TEST: Example System"):
print("Expected solutions:"):
print("  x1 = -(y1^2-2*y1+3)/(y1*y2-1)"):
print("  x2 = -(y1*y2-y1-3*y2+1)/(y1*y2-1)"):
print(""):

# Sys, Vars, params, p := get_data("example"):
# Sys, Vars, params, p := get_data("small_Sys"):
Sys, Vars, params, p := get_data("small_sys_low_deg"):
# Sys, Vars, params, p := get_data("bspline"):

print("System:", Sys):
print("Variables:", Vars):
print("Parameters:", params):
print("Prime:", p):

# Create black box
LBB := Constuct_Sys_Blackbox(Sys, Vars, params):

# Test black box
test_pt := [2, 3]:
print(""):
print("Testing black box at point", test_pt):
result := LBB(test_pt, p):
print("Result:", result):

# Expected values
den_val := (2*3-1) mod p:  # = 5
num1_val := -(2^2-2*2+3) mod p:  # = -(4-4+3) = -3 = 104 mod 107
num2_val := -(2*3-2-3*3+1) mod p:  # = -(6-2-9+1) = -(-4) = 4
exp1 := num1_val/den_val mod p:  # 104/5 mod 107
exp2 := num2_val/den_val mod p:   # 4/5 mod 107
print("Expected: [", exp1, ",", exp2, "]"):

# Run MRFI
print(""):
print("Running MRFI_Vector..."):
print(""):

num_vars := nops(params):
num_eqn := nops(Vars):

try
    Numerators, Denominators := MRFI_Vector(LBB, num_vars, num_eqn, params, p):
    
    print(""):
    print("========================================"):
    print("RESULTS:"):
    print("========================================"):
    
    for k from 1 to num_eqn do
        print(""):
        print("Component", k, ":");
        print("  Numerator:", Numerators[k]):
        print("  Denominator:", Denominators[k]):
    end do:
    
    # Expected polynomials
    print(""):
    print("EXPECTED:"):
    expected_num1 := -(y1^2-2*y1+3) mod p:
    expected_num2 := -(y1*y2-y1-3*y2+1) mod p:
    expected_den := y1*y2-1 mod p:
    print("  Numerator 1:", expected_num1):
    print("  Numerator 2:", expected_num2):
    print("  Denominator:", expected_den):
    
    # Verification
    print(""):
    print("VERIFICATION:"):
    print("-------------"):
    
    r := rand(p):
    num_tests := 5:
    success_count := 0:
    
    for test_num from 1 to num_tests do
        test_point := [r(), r()]:
        
        # Skip if denominator is zero
        den_check := Eval(expected_den, {y1=test_point[1], y2=test_point[2]}) mod p:
        if den_check = 0 then
            print("Skipping test", test_num, "(denominator zero)"):
            next:
        end if:
        
        print(""):
        print("Test", test_num, "at point", test_point):
        
        # Evaluate black box
        bb_vals := LBB(test_point, p):
        print("  Black box:", bb_vals):
        
        # Evaluate recovered functions
        rec_vals := []:
        for k from 1 to num_eqn do
            n_eval := Eval(Numerators[k], {y1=test_point[1], y2=test_point[2]}) mod p:
            d_eval := Eval(Denominators[k], {y1=test_point[1], y2=test_point[2]}) mod p:
            if d_eval <> 0 then
                rec_vals := [op(rec_vals), n_eval/d_eval mod p]:
            else
                print("  Error: Recovered denominator zero!"):
                break:
            end if:
        end do:
        
        if nops(rec_vals) = num_eqn then
            print("  Recovered:", rec_vals):
            
            # Check match
             differ := [seq(bb_vals[i] - rec_vals[i] mod p, i=1..num_eqn)]:
            if convert(differ, set) = {0} then
                print("  ✓ MATCH"):
                success_count := success_count + 1:
            else
                print("  ✗ MISMATCH, differ =", differ):
            end if:
        end if:
    end do:
    
    print(""):
    print("========================================"):
    print("Success rate:", success_count, "/", num_tests):
    print("========================================"):
    
catch:
    print("ERROR:", lasterror()):
end try:
