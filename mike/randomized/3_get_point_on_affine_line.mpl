
########################################################
# 3. GET POINT ON AFFINE LINE PSI
########################################################

get_point_on_affine_line := proc(num_var::posint, alpha::list, beta_::list,sigma_::list, p::prime, T::posint)
    description "Generates evaluation points on a parametric affine line in finite field Z_p^n":
    
    option remember;  # Cache results for repeated calls with same parameters
    
# 1. Parameters:
# -----------
# num_var : positive integer
#     The dimension of the affine space (number of variables in the polynomial system)
#     Constraint: num_var >= 2
#
# alpha : list of integers. These are points in Z_p. Map Psi will find the points on the line for these alphas.
#     List of T parameter values for the affine line, typically [alpha_1, alpha_2, ..., alpha_T]
#     These serve as the univariate parameter for the line parameterization
#     Requirement: nops(alpha) >= T, all elements in range [0, p-1]
#
# beta_ : list of integers
#     Direction vector coefficients [beta_1, beta_2, ..., beta_{num_var-1}]
#     Defines the direction of the affine line in the finite field space
#     Requirement: nops(beta_) = num_var - 1, all elements in range [0, p-1]
#
# sigma_ : list of integers
#     Base point coordinates [sigma_1, sigma_2, ..., sigma_{num_var}]
#     Defines the initial point through which the affine line passes
#     Requirement: nops(sigma_) = num_var, all elements in range [0, p-1]
#
# p : prime number
#     The characteristic of the finite field Z_p
#     All arithmetic operations are performed modulo p
#
# T : positive integer
#     Number of evaluation points to generate on the affine line
#     Constraint: T <= p (to ensure distinct points)
#
# 2. Returns:
# --------
# list of lists
#     A list of T points, where each point on the affine line is a list
#     Format: [[psi_1,1, psi_1,2, ..., psi_1,n],...,[psi_T,1, psi_T,2, ..., psi_T,n]]
#     where psi_i,j represents the j-th coordinate of the i-th point
#
# 3. Mathematical Description:
# ------------------------
# This procedure computes points on an affine line in Z_p^n parameterized as:
#     psi(alpha) = sigma + alpha * beta
# where:
#     - sigma = (sigma_1, sigma_2, ..., sigma_n) is the base point
#     - beta is the direction vector
#     - alpha is the line parameter
#
# Specifically, for each parameter value alpha_i:
#     psi_i,1 = alpha_i (first coordinate uses the parameter directly)
#     psi_i,j = beta_{j-1} * (alpha_i - sigma_1) + sigma_j (mod p) for j = 2, ..., num_var
#
# This creates a line that:
#     1. Passes through point sigma when alpha = sigma_1
#     2. Has direction determined by the beta coefficients
#     3. Ensures the first coordinate equals the parameter value
#
# 4. Algorithm Complexity:
# --------------------
# Time Complexity: O(T * num_var)
# Space Complexity: O(T * num_var)
#
# 5. Error Conditions:
# ----------------
# - If nops(alpha) < T: May cause index out of bounds error
# - If nops(beta_) != num_var - 1: Incorrect parameterization
# - If nops(sigma_) != num_var: Incorrect base point specification
#
# 6. Notes:
# ------
# The affine line construction ensures good separation of evaluation points
#
# 7. References:
# ----------

    # [1] Kaltofen, E. & Yang, Z. (2007). "Sparse Multivariate Function Recovery"
    
    local psi, nv, np, i;
    global num_lines:
    num_lines:=num_lines+1:
    # Checking preconditions
    if num_var < 2 then
        error "Number of variables must be at least 2, got %1", num_var;
    end if;
    if nops(alpha) < T then
        error "Insufficient alpha values: need %1, got %2", T, nops(alpha);
    end if;
    # Input validation (recommended to add)
    ASSERT(num_var >= 2, "Number of variables must be at least 2");
    ASSERT(nops(alpha) >= T, "Insufficient alpha values for T points");
    ASSERT(nops(beta_) = num_var - 1, "beta_ must have exactly num_var-1 elements");
    ASSERT(nops(sigma_) >= num_var, "sigma_ must have at least num_var elements");
    
    
    # Generate T points on the affine line
    for np from 1 to T do 
        # First coordinate is always the parameter value
        psi[np][1] := alpha[np];
        
        # Remaining coordinates follow the affine line formula
        for nv from 2 to num_var do 
            psi[np][nv] := beta_[nv-1]*alpha[np] - beta_[nv-1]*sigma_[1] + sigma_[nv] mod p;
        end do;
    end do;
    
    # Convert the array of points to a list of lists
    return [seq(convert(psi[i], list), i=1..T)];
end proc: