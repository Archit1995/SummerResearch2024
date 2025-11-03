
MRFI_Vector:=proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, Primes, shift_, sigma_, num_list, den_list, 
          numerator_done, denominator_done, T, T_old,
          mqrfr_results, lin_sys, num_points_mqrfr,
          F_current, G_current, deg_num, deg_den,
          lambda_num_list, lambda_den_list, terms_num_list, terms_den_list,
          R_num_list, R_den_list, Roots_num_list, Roots_den_list,
          num_mono_list, den_mono_list, coeff_num_list, coeff_den_list,
          final_num_list, final_den_list, r_, temp, common_den_flag,
          den_lc, den_lc_inv:
    
    lprint("MRFI ========================================"):
    lprint("MRFI Starting MRFI_Vector"):
    lprint("MRFI Number of parameters:", num_vars):
    lprint("MRFI Number of equations:", num_eqn):
    lprint("MRFI ========================================"):
    
    r_ := rand(p):
    Primes := [seq(ithprime(i), i=1..num_vars)]:
    shift_ := [seq(r_(), i=1..num_vars-1)]:
    sigma_ := []:
    
    # Initialize storage
    num_list := [seq([], i=1..num_eqn)]:
    den_list := [seq([], i=1..num_eqn)]:
    numerator_done := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    lprint("MRFI numerator_done: ",numerator_done):
    lprint("MRFI denominator_done: ",denominator_done):
    T := 4:
    T_old := 1:
    
    # Initial NDSA call
    mqrfr_results, lin_sys := NDSA(B, [seq(1, i=1..num_vars)], shift_, num_vars, p, T):
        lprint("MRFI numerator_done: ",numerator_done):
        lprint("MRFI denominator_done: ",denominator_done):
    if not lin_sys then
        error "Expected a linear system (vector output)":
    end if:
    
    # Process initial results
    F_current := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    lprint("MRFI F_current: ",F_current):
    G_current := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
    lprint("MRFI G_current: ",G_current):

    for k from 1 to num_eqn do
        num_list[k] := [eval(F_current[k], x=1)]:
        
        den_list[k] := [eval(G_current[k], x=1)]:
    end do:
    lprint("MRFI num_list: ",num_list):
    lprint("MRFI den_list: ",den_list):
    deg_num := [seq(degree(F_current[i], x), i=1..nops(F_current))]:
    lprint("MRFI deg_num: ",deg_num):
    deg_den := [seq(degree(G_current[i], x), i=1..nops(G_current))]:
    lprint("MRFI deg_den: ",deg_den):
    #     num_points_mqrfr:=deg_num+deg_den+2:
    num_points_mqrfr := max(op(deg_num)) + max(op(deg_den)) + 2:
    lprint("MRFI num_points_mqrfr: ",num_points_mqrfr):
    # Check for common denominator
    common_den_flag := true:
    lprint("MRFI G_current: ",G_current):
    for k from 2 to num_eqn do
        if G_current[k] <> G_current[1] then
            common_den_flag := false:
            break:
        end if:
    end do:
    lprint("MRFI Common denominator:", common_den_flag):
    print("______________________________________________________________________________"):
    print("in main evaluation loop of MRFI"):
    
    # Main evaluation loop
    while true do
        lprint("MRFI T_old=", T_old):
        lprint("MRFI T=", T):
        
        for j from T_old to 2*T-1 do
            sigma_ := [op(sigma_), [seq(Primes[i]^j mod p, i=1..nops(Primes))]]:
            
            # mqrfr_results, temp, lin_sys := NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr):
            # Add an index of remaining here
            mqrfr_results, lin_sys := NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr):
            
            F_current := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
            lprint("MRFI F_current: ",F_current):
            G_current := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
            lprint("MRFI G_current: ",G_current):
            print("______________________________________________________________________________"):
            print("Udating num_list and den_list at sigma_[",j,"]=",sigma_[j]):


            for k from 1 to num_eqn do
                lprint("MRFI k=",k):
                lprint("MRFI num_eqn=",num_eqn):
                lprint("MRFI numerator_done[",k,"]=",numerator_done[k]):
                lprint("MRFI num_list: ",num_list):
                lprint("MRFI den_list: ",den_list):
                if not numerator_done[k] then
                    num_list[k] := [op(num_list[k]), eval(F_current[k], x=sigma_[j][1]) mod p]:
                end if:
                if not denominator_done[k] then
                    den_list[k] := [op(den_list[k]), eval(G_current[k], x=sigma_[j][1]) mod p]:
                end if:
            end do:
            print("MRFI num_list after loop: ",num_list):
        end do:
        lprint("MRFI num_list: ",num_list):
        
        # Apply BMEA
        lambda_num_list := []:
        terms_num_list := []:
        R_num_list := []:
        
        all_num_done := true:
        for k from 1 to num_eqn do
            lprint("MRFI k=",k):
            lprint("MRFI num_eqn=",num_eqn):
            lprint("MRFI numerator_done[",k,"]=",numerator_done[k]): 
            numerator_done[k] := false: 
            if not numerator_done[k] then
                temp := BMEA(num_list[k], p, Z):
                lprint("MRFI temp (numerator)=",temp):
                lambda_num_list := [op(lambda_num_list), temp]:
                lprint("MRFI lambda_num_list: ",lambda_num_list):
                terms_num_list := [op(terms_num_list), degree(temp, Z)]:
                lprint("MRFI terms_num_list: ",terms_num_list):
                R_num_list := [op(R_num_list), Roots(temp) mod p]:
                lprint("MRFI R_num_list: ",R_num_list):
                print("--------------------------------------------------------------------------------"):
                print("Checking termination condition for numerator"):


                
                print("nops(R_num_list[",k,"]): ", nops(R_num_list[k])):
                print("R_num_list[",k,"][1][1]: ", R_num_list[k][1][1]);
                print("terms_num_list[",k,"]: ", terms_num_list[k]);
                print("T: ", T);



                if nops(R_num_list[k]) > 0 and R_num_list[k][1][1] = 0 then
                    R_num_list[k] := remove(x->x=[0,1], R_num_list[k]):
                    terms_num_list[k] := terms_num_list[k] - 1:
                end if:
                # if nops(R_num_list[k]) = terms_num_list[k] and terms_num_list[k] < T then
                if nops(R_num_list[k]) = terms_num_list[k] and terms_num_list[k] <= T then
                    numerator_done[k] := true:
                else
                    all_num_done := false:
                end if:
            end if:
        end do:
        print("===================================================================================="):
        print("MFRI processing denominators now"):
        # Process denominators
        lprint("MRFI den_list: ",den_list):
        lambda_den_list := []:
        terms_den_list := []:
        R_den_list := []:
        
        all_den_done := true:
        lprint("MRFI common_den_flag: ",common_den_flag):
        if common_den_flag then
            lprint("MRFI denominator_done: ",denominator_done):
            denominator_done[1] := false:
            if not denominator_done[1] then
                temp := BMEA(den_list[1], p, Z):
                lprint("MRFI temp (denominator)=",temp):
                lambda_den_list := [temp]:
                lprint("MRFI lambda_den_list: ",lambda_den_list):
                terms_den_list := [degree(temp, Z)]:
                lprint("MRFI terms_den_list: ",terms_den_list):
                R_den_list := [Roots(temp) mod p]:
                lprint("MRFI R_den_list: ",R_den_list):
                
                if nops(R_den_list[1]) > 0 and R_den_list[1][1][1] = 0 then
                    R_den_list[1] := remove(x->x=[0,1], R_den_list[1]):
                    terms_den_list[1] := terms_den_list[1] - 1:
                end if:
                
                if nops(R_den_list[1]) = terms_den_list[1] and terms_den_list[1] < T then
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do:
                else
                    all_den_done := false:
                end if:
            end if:
        else
            for k from 1 to num_eqn do
                if not denominator_done[k] then
                    temp := BMEA(den_list[k], p, Z):
                    lambda_den_list := [op(lambda_den_list), temp]:
                    terms_den_list := [op(terms_den_list), degree(temp, Z)]:
                    R_den_list := [op(R_den_list), Roots(temp) mod p]:
                    print("--------------------------------------------------------------------------------"):
                    print("Checking termination condition for denominator"):
                    print("nops(R_den_list[",k,"]): ", nops(R_den_list[k])):
                    print("R_den_list[",k,"][1][1]: ", R_den_list[k][1][1]);
                    print("terms_den_list[",k,"]: ", terms_den_list[k]);
                    print("T: ", T);
                    
                    if nops(R_den_list[k]) > 0 and R_den_list[k][1][1] = 0 then
                        R_den_list[k] := remove(x->x=[0,1], R_den_list[k]):
                        terms_den_list[k] := terms_den_list[k] - 1:
                    end if:
                    
                    if nops(R_den_list[k]) = terms_den_list[k] and terms_den_list[k] < T then
                        denominator_done[k] := true:
                    else
                        all_den_done := false:
                    end if:
                end if:
            end do:
        end if:
        
        if all_num_done and all_den_done then
            lprint("MRFI All components recovered!"):
            break:
        end if:
        
        T_old := 2*T:
        T := T*2:
        print("______________________________________________________________________________"):
    end do:
    
    # Extract roots and build polynomials
    Roots_num_list := [seq([seq(r[1], r in R_num_list[k])], k=1..num_eqn)]:
    
    if common_den_flag then
        Roots_den_list := [[seq(r[1], r in R_den_list[1])]]:
        for k from 2 to num_eqn do
            Roots_den_list := [op(Roots_den_list), Roots_den_list[1]]:
        end do:
    else
        Roots_den_list := [seq([seq(r[1], r in R_den_list[k])], k=1..num_eqn)]:
    end if:
    
    # Generate monomials and coefficients
    num_mono_list := []:
    coeff_num_list := []:
    final_num_list := []:
    
    for k from 1 to num_eqn do
        temp := generate_monomials(Roots_num_list[k], num_vars, Primes, vars):
        if temp = FAIL then return FAIL: end if:
        num_mono_list := [op(num_mono_list), temp]:
        
        coeff_num_list := [op(coeff_num_list), 
            Zippel_Vandermonde_solver(num_list[k], terms_num_list[k], 
                                     Roots_num_list[k], lambda_num_list[k], p)]:
        
        final_num_list := [op(final_num_list), 
            construct_final_polynomial(coeff_num_list[k], num_mono_list[k])]:
    end do:
    
    # Generate denominators
    den_mono_list := []:
    coeff_den_list := []:
    final_den_list := []:
    
    if common_den_flag then
        temp := generate_monomials(Roots_den_list[1], num_vars, Primes, vars):
        if temp = FAIL then return FAIL: end if:
        
        coeff_den := Zippel_Vandermonde_solver(den_list[1], terms_den_list[1],
                                               Roots_den_list[1], lambda_den_list[1], p):
        
        temp_den := construct_final_polynomial(coeff_den, temp):
        final_den_list := [seq(temp_den, k=1..num_eqn)]:
    else
        for k from 1 to num_eqn do
            temp := generate_monomials(Roots_den_list[k], num_vars, Primes, vars):
            if temp = FAIL then return FAIL: end if:
            
            coeff_den := Zippel_Vandermonde_solver(den_list[k], terms_den_list[k],
                                                   Roots_den_list[k], lambda_den_list[k], p):
            
            final_den_list := [op(final_den_list), 
                construct_final_polynomial(coeff_den, temp)]:
        end do:
    end if:
    
    # Normalize
    for k from 1 to num_eqn do
        den_lc := lcoeff(final_den_list[k], order=grlex(seq(vars[i], i=1..num_vars))):
        den_lc_inv := 1/den_lc mod p:
        final_num_list[k] := final_num_list[k] * den_lc_inv mod p:
        final_den_list[k] := final_den_list[k] * den_lc_inv mod p:
    end do:
    
    lprint("MRFI ========================================"):
    lprint("MRFI RECOVERY COMPLETE"):
    lprint("MRFI ========================================"):
    
    return final_num_list, final_den_list:
end proc:
