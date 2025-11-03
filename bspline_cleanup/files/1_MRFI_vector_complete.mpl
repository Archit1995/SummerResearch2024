MRFI_Vector:=proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, Primes, shift_, sigma_, num_list, den_list, 
          numerator_done, denominator_done, T, T_old,
          mqrfr_results, lin_sys, num_points_mqrfr,
          F_current, G_current, deg_num, deg_den,
          lambda_num_list, lambda_den_list, terms_num_list, terms_den_list,
          R_num_list, R_den_list, Roots_num_list, Roots_den_list,
          num_mono_list, den_mono_list, coeff_num_list, coeff_den_list,
          final_num_list, final_den_list, r_, temp, common_den_flag;
    
    print("In MRFI_Vector"):
    print("Number of equations:", num_eqn):
    print("Number of variables:", num_vars):
    
    r_ := rand(p):
    Primes := [seq(ithprime(i), i=1..num_vars)]:
    shift_ := [seq(r_(), i=1..num_vars-1)]:
    print("shift =", shift_):
    
    sigma_ := []:
    
    # Initialize storage for evaluations of each component
    # num_list[i] stores evaluations of numerator i
    num_list := [seq([], i=1..num_eqn)]:
    den_list := [seq([], i=1..num_eqn)]:
    
    # Track completion status for each component
    numerator_done := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    
    T := 4:
    T_old := 1:
    
    # Initial call to NDSA to determine degrees
    print("Initial NDSA call"):
    mqrfr_results, num_points_mqrfr, lin_sys := NDSA(B, [seq(1, i=1..num_vars)], shift_, num_vars, p, T):
    
    if not lin_sys then
        error "Expected a linear system (vector output), but got scalar":
    end if:
    
    print("Number of components:", nops(mqrfr_results)):
    
    # Process initial results
    F_current := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    G_current := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
    
    # Store initial evaluations
    for k from 1 to num_eqn do
        num_list[k] := [op(num_list[k]), eval(F_current[k], x=1)]:
        den_list[k] := [op(den_list[k]), eval(G_current[k], x=1)]:
    end do:
    
    # Get degrees
    deg_num := [seq(degree(F_current[i], x), i=1..nops(F_current))]:
    deg_den := [seq(degree(G_current[i], x), i=1..nops(G_current))]:
    print("Numerator degrees:", deg_num):
    print("Denominator degrees:", deg_den):
    
    # Calculate required points
    num_points_mqrfr := max(op(deg_num)) + max(op(deg_den)) + 2:
    print("num_points_mqrfr =", num_points_mqrfr):
    
    # Check if all denominators are the same (common denominator case)
    common_den_flag := true:
    for k from 2 to num_eqn do
        if degree(G_current[k], x) <> degree(G_current[1], x) then
            common_den_flag := false:
            break:
        end if:
    end do:
    print("Common denominator detected:", common_den_flag):
    
    # Main loop for gathering evaluations
    while true do
        print("T=", T, "T_old=", T_old):
        
        for j from T_old to 2*T-1 do
            sigma_ := [op(sigma_), [seq(Primes[i]^j mod p, i=1..nops(Primes))]]:
            print("sigma_[", j, "]=", sigma_[j]):
            
            # Get evaluations at this point
            mqrfr_results, temp, lin_sys := NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr):
            
            if not lin_sys then
                error "Expected vector output from NDSA":
            end if:
            
            # Extract current F and G
            F_current := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
            G_current := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
            
            # Evaluate and store for each component
            for k from 1 to num_eqn do
                if not numerator_done[k] then
                    num_list[k] := [op(num_list[k]), eval(F_current[k], x=sigma_[j][1]) mod p]:
                end if:
                if not denominator_done[k] then
                    den_list[k] := [op(den_list[k]), eval(G_current[k], x=sigma_[j][1]) mod p]:
                end if:
            end do:
        end do:
        
        # Apply BMEA to each component
        print("Applying BMEA..."):
        
        # Process numerators
        lambda_num_list := []:
        terms_num_list := []:
        R_num_list := []:
        
        for k from 1 to num_eqn do
            if not numerator_done[k] then
                print("Processing numerator", k):
                lambda_num_list := [op(lambda_num_list), BMEA(num_list[k], p, Z)]:
                terms_num_list := [op(terms_num_list), degree(lambda_num_list[k], Z)]:
                R_num_list := [op(R_num_list), Roots(lambda_num_list[k]) mod p]:
                
                # Remove zero root if present
                if nops(R_num_list[k]) > 0 and R_num_list[k][1][1] = 0 then
                    R_num_list[k] := remove(x->x=[0,1], R_num_list[k]):
                    terms_num_list[k] := terms_num_list[k] - 1:
                    lambda_num_list[k] := quo(lambda_num_list[k], Z, Z) mod p:
                end if:
                
                # Check success
                if nops(R_num_list[k]) = terms_num_list[k] and terms_num_list[k] < T then
                    print("Numerator", k, "Success!"):
                    numerator_done[k] := true:
                end if:
            end if:
        end do:
        
        # Process denominators
        lambda_den_list := []:
        terms_den_list := []:
        R_den_list := []:
        
        if common_den_flag then
            # Use only first denominator for common case
            if not denominator_done[1] then
                print("Processing common denominator"):
                lambda_den_list := [BMEA(den_list[1], p, Z)]:
                terms_den_list := [degree(lambda_den_list[1], Z)]:
                R_den_list := [Roots(lambda_den_list[1]) mod p]:
                
                # Remove zero root if present
                if nops(R_den_list[1]) > 0 and R_den_list[1][1][1] = 0 then
                    R_den_list[1] := remove(x->x=[0,1], R_den_list[1]):
                    terms_den_list[1] := terms_den_list[1] - 1:
                    lambda_den_list[1] := quo(lambda_den_list[1], Z, Z) mod p:
                end if:
                
                # Check success
                if nops(R_den_list[1]) = terms_den_list[1] and terms_den_list[1] < T then
                    print("Common denominator Success!"):
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do:
                end if:
            end if:
        else
            # Process each denominator separately
            for k from 1 to num_eqn do
                if not denominator_done[k] then
                    print("Processing denominator", k):
                    lambda_den_list := [op(lambda_den_list), BMEA(den_list[k], p, Z)]:
                    terms_den_list := [op(terms_den_list), degree(lambda_den_list[k], Z)]:
                    R_den_list := [op(R_den_list), Roots(lambda_den_list[k]) mod p]:
                    
                    # Remove zero root if present
                    if nops(R_den_list[k]) > 0 and R_den_list[k][1][1] = 0 then
                        R_den_list[k] := remove(x->x=[0,1], R_den_list[k]):
                        terms_den_list[k] := terms_den_list[k] - 1:
                        lambda_den_list[k] := quo(lambda_den_list[k], Z, Z) mod p:
                    end if:
                    
                    # Check success
                    if nops(R_den_list[k]) = terms_den_list[k] and terms_den_list[k] < T then
                        print("Denominator", k, "Success!"):
                        denominator_done[k] := true:
                    end if:
                end if:
            end do:
        end if:
        
        # Check if all components are done
        if andmap(x->x, numerator_done) and andmap(x->x, denominator_done) then
            print("All components successfully recovered!"):
            break:
        end if:
        
        # Double the bound for next iteration
        T_old := 2*T:
        T := T*2:
    end do:
    
    # Extract roots
    Roots_num_list := []:
    for k from 1 to num_eqn do
        Roots_num_list := [op(Roots_num_list), [seq(r[1], r in R_num_list[k])]]:
    end do:
    
    Roots_den_list := []:
    if common_den_flag then
        Roots_den_list := [[seq(r[1], r in R_den_list[1])]]:
        # Replicate for all components
        for k from 2 to num_eqn do
            Roots_den_list := [op(Roots_den_list), Roots_den_list[1]]:
        end do:
    else
        for k from 1 to num_eqn do
            Roots_den_list := [op(Roots_den_list), [seq(r[1], r in R_den_list[k])]]:
        end do:
    end if:
    
    print("---------------------------------------------"):
    for k from 1 to num_eqn do
        print("Roots_num[", k, "]:", Roots_num_list[k]):
    end do:
    if common_den_flag then
        print("Common Roots_den:", Roots_den_list[1]):
    else
        for k from 1 to num_eqn do
            print("Roots_den[", k, "]:", Roots_den_list[k]):
        end do:
    end if:
    print("---------------------------------------------"):
    
    # Generate monomials for each numerator
    num_mono_list := []:
    for k from 1 to num_eqn do
        temp := generate_monomials(Roots_num_list[k], num_vars, Primes, vars):
        if temp = FAIL then
            return FAIL:
        end if:
        num_mono_list := [op(num_mono_list), temp]:
        print("num_mono[", k, "]:", temp):
    end do:
    
    # Generate monomials for denominators
    den_mono_list := []:
    if common_den_flag then
        temp := generate_monomials(Roots_den_list[1], num_vars, Primes, vars):
        if temp = FAIL then
            return FAIL:
        end if:
        den_mono_list := [temp]:
        for k from 2 to num_eqn do
            den_mono_list := [op(den_mono_list), temp]:
        end do:
    else
        for k from 1 to num_eqn do
            temp := generate_monomials(Roots_den_list[k], num_vars, Primes, vars):
            if temp = FAIL then
                return FAIL:
            end if:
            den_mono_list := [op(den_mono_list), temp]:
        end do:
    end if:
    
    # Compute coefficients using Vandermonde solver
    coeff_num_list := []:
    for k from 1 to num_eqn do
        coeff_num_list := [op(coeff_num_list), 
            Zippel_Vandermonde_solver(num_list[k], terms_num_list[k], 
                                     Roots_num_list[k], lambda_num_list[k], p)]:
        print("coeff_num[", k, "]:", coeff_num_list[k]):
    end do:
    
    coeff_den_list := []:
    if common_den_flag then
        coeff_den_list := [Zippel_Vandermonde_solver(den_list[1], terms_den_list[1],
                                                     Roots_den_list[1], lambda_den_list[1], p)]:
        for k from 2 to num_eqn do
            coeff_den_list := [op(coeff_den_list), coeff_den_list[1]]:
        end do:
    else
        for k from 1 to num_eqn do
            coeff_den_list := [op(coeff_den_list),
                Zippel_Vandermonde_solver(den_list[k], terms_den_list[k],
                                         Roots_den_list[k], lambda_den_list[k], p)]:
        end do:
    end if:
    
    # Construct final polynomials
    final_num_list := []:
    for k from 1 to num_eqn do
        final_num_list := [op(final_num_list), 
            construct_final_polynomial(coeff_num_list[k], num_mono_list[k])]:
    end do:
    
    final_den_list := []:
    if common_den_flag then
        temp := construct_final_polynomial(coeff_den_list[1], den_mono_list[1]):
        final_den_list := [seq(temp, k=1..num_eqn)]:
    else
        for k from 1 to num_eqn do
            final_den_list := [op(final_den_list),
                construct_final_polynomial(coeff_den_list[k], den_mono_list[k])]:
        end do:
    end if:
    
    # Normalize by leading coefficients
    print("==============================================="):
    print("FINAL RESULTS:"):
    print("==============================================="):
    
    for k from 1 to num_eqn do
        den_lc := lcoeff(final_den_list[k], order=grlex(seq(vars[i], i=1..num_vars))):
        den_lc_inv := 1/den_lc mod p:
        
        final_num_list[k] := final_num_list[k] * den_lc_inv mod p:
        final_den_list[k] := final_den_list[k] * den_lc_inv mod p:
        
        print("Component", k, ":"):
        print("  Numerator:", final_num_list[k]):
        print("  Denominator:", final_den_list[k]):
    end do:
    
    return final_num_list, final_den_list:
end proc:
