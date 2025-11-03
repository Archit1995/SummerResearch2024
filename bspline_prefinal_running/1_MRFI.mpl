
MRFI_Vector:=proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, Primes, shift_, sigma_, num_eval, den_eval, 
          numerator_done, denominator_done, T, T_old,
          mqrfr_results, lin_sys, num_points_mqrfr,
          Numerators, Denominiators, deg_num, deg_den,
          lambda_num_eval, lambda_den_eval, terms_num_eval, terms_den_eval,
          R_num_eval, R_den_eval, Roots_num_eval, Roots_den_eval,
          num_mono_list, den_mono_list, coeff_num_eval, coeff_den_eval,
          final_num_eval, final_den_eval, r_, temp, common_den_flag,
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
    
    # Initializing accumulators and flags
    num_eval := [seq([], i=1..num_eqn)]:
    den_eval := [seq([], i=1..num_eqn)]:

    # lambda_num_eval := [seq([], i=1..num_eqn)]:
    # terms_num_eval := [seq([], i=1..num_eqn)]:
    # R_num_eval := [seq([], i=1..num_eqn)]:
    # lambda_den_eval := [seq([], i=1..num_eqn)]:
    # terms_den_eval := [seq([], i=1..num_eqn)]:
    # R_den_eval := [seq([], i=1..num_eqn)]:
    
    numerator_done := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    bmea_done := [seq(false, i=1..num_eqn)]:
    all_done := true:

    # lambda_num_eval := []:
    # terms_num_eval := []:
    # R_num_eval := []:
    # lambda_den_eval := []:
    # terms_den_eval := []:
    # R_den_eval := []:
    
    lambda_num_eval := table():
    terms_num_eval := table():
    R_num_eval := table():
    lambda_den_eval := table():
    terms_den_eval := table():
    R_den_eval := table():

    # Initialize all entries from 1 to 21
    for i from 1 to num_eqn do
        lambda_num_eval[i] := []:
        terms_num_eval[i] := []:
        R_num_eval[i] := []:
        lambda_den_eval[i] := []:
        terms_den_eval[i] := []:
        R_den_eval[i] := []:
    end do:

    lprint("MRFI numerator_done: ",numerator_done):
    lprint("MRFI denominator_done: ",denominator_done):
    T := 4:
    T_old := 1:
    num_points_mqrfr:=[seq(0, i=1..num_eqn)]:
    
    # Initial NDSA call
    # mqrfr_results, num_points_mqrfr, lin_sys := NDSA(B, [seq(1, i=1..num_vars)], shift_, num_vars, p, T,num_eqn):
    #mqrfr_results = [N(x),D(x),deg(qmax),lc(D(x))]
    mqrfr_results, lin_sys := NDSA(B, [seq(1, i=1..num_vars)], shift_, num_vars, p, T,num_eqn):
        # lprint("MRFI denominator_done: ",denominator_done):
    
    if not lin_sys then
        error "Expected a linear system (vector output)":
    end if:
    
    # Process initial results
    Numerators := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    # lprint("MRFI Numerators: ",Numerators):
    Denominiators := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
    # lprint("MRFI Denominiators: ",Denominiators):
    # print("num_eqn: ", num_eqn):
    for k from 1 to num_eqn do
        # print("Evaluating univarite numerators and denominators at sigma[1] and storing them in accumulators ", k):
        num_eval[k] := [eval(Numerators[k], x=1)]:
        # print("num_eval[",k,"]=",num_eval[k]):
        den_eval[k] := [eval(Denominiators[k], x=1)]:
        # print("den_eval[",k,"]=",den_eval[k]):
    end do:
    # print("MRFI num_eval:=[[N1(1),N1(2),N1(4)...],[N2(1),N2(2),N2(4)...]] ",num_eval):
    # print("MRFI den_eval:=[[D1(1),D1(2),D1(4)...],[D2(1),D2(2),D2(4)...]] ",den_eval):
    deg_num := [seq(degree(Numerators[i], x), i=1..nops(Numerators))]:
    # lprint("MRFI deg_num: ",deg_num):
    deg_den := [seq(degree(Denominiators[i], x), i=1..nops(Denominiators))]:
    # lprint("MRFI deg_den: ",deg_den):
    #     num_points_mqrfr:=deg_num+deg_den+2:

    for i from 1 to numelems(deg_den) do 
        num_points_mqrfr[i] := deg_num[i] + deg_den[i] + 2:
    end do:
    max_num_points_mqrfr := max(op(deg_num)) + max(op(deg_den)) + 2:
    print("MRFI num_points_mqrfr: ",num_points_mqrfr):
    print("MRFI max_num_points_mqrfr: ",max_num_points_mqrfr):
    # Check for common denominator
    common_den_flag := true:
   
    for k from 2 to num_eqn do
        if Denominiators[k] <> Denominiators[1] then
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
            # print("MRFI Evaluating at sigma_[",nops(sigma_),"]=",sigma_[nops(sigma_)]):
            # mqrfr_results, temp, lin_sys := NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr):
            # Add an index of remaining here

            mqrfr_results:= Deterministic_NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr,max_num_points_mqrfr,num_eqn):

            # mqrfr_results, lin_sys := NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr,num_eqn):
            # print("MRFI NDSA results at sigma_[",j,"]=",sigma_[j],": ",mqrfr_results):
            Numerators := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
            # lprint("MRFI Numerators: ",Numerators):
            Denominiators := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
            # lprint("MRFI Denominiators: ",Denominiators):
            print("______________________________________________________________________________"):
            # print("Updating num_eval and den_eval at sigma_[",j,"]=",sigma_[j]):


            for k from 1 to num_eqn do
                # print("MRFI k=",k):
                # print("MRFI num_eqn=",num_eqn):
                # print("MRFI numerator_done[",k,"]=",numerator_done[k]):
                # print("MRFI num_eval: ",num_eval):
                # print("MRFI den_eval: ",den_eval):
                if not numerator_done[k] then
                    num_eval[k] := [op(num_eval[k]), eval(Numerators[k], x=sigma_[j][1]) mod p]:
                end if:
                if not denominator_done[k] then
                    den_eval[k] := [op(den_eval[k]), eval(Denominiators[k], x=sigma_[j][1]) mod p]:
                end if:
                # print("--------------------------------------------------------------------------------"):
            end do:
            # print("MRFI num_eval after loop: ",num_eval):
            # print("MRFI den_eval after loop: ",den_eval):
            # print("______________________________________________________________________________"):
        end do:
        print("MRFI num_eval: ",num_eval):
        print("MRFI den_eval: ",den_eval):
        print("______________________________________________________________________________"):
        print("numerator_done: ",numerator_done):
        print("denominator_done: ",denominator_done):
        print("--------------------------------------------------------------------------------"):
        # break:
        # Apply BMEA
        lambda_num_eval := []:
        terms_num_eval := []:
        R_num_eval := []:
        # lambda_den_eval := []:
        # terms_den_eval := []:
        # R_den_eval := []:
        all_done := true:
        for k from 1 to num_eqn do
            lprint("MRFI k=",k):
            lprint("MRFI numerator_done[",k,"]=",numerator_done[k]):
            if  numerator_done[k] then print("MRFI Skipping BMEA for numerator component ",k," as already done"): next: end if;
                # temp := BMEA(num_eval[k], p, Z):
                # lprint("MRFI temp (numerator)=",temp):
                # accumulator list 
                
                # lambda_num_eval:= [op(lambda_num_eval), temp]:
                # terms_num_eval:= [op(terms_num_eval), degree(temp, Z)]:
                # R_num_eval := [op(R_num_eval), Roots(temp) mod p]:

                lambda_num_eval := [op(lambda_num_eval), BMEA(num_eval[k], p, Z)]:
                terms_num_eval:= [op(terms_num_eval), degree(lambda_num_eval[k], Z)]:
                R_num_eval := [op(R_num_eval), Roots(lambda_num_eval[k]) mod p]:
                
                
                # hash
                # lambda_num_eval[k] := [temp]:
                # terms_num_eval[k] := degree(temp, Z):
                # R_num_eval[k] := Roots(temp) mod p:

                lprint("MRFI lambda_num_eval: ",lambda_num_eval[k]):
                lprint("MRFI terms_num_eval: ",terms_num_eval[k]):
                lprint("MRFI R_num_eval: ",R_num_eval[k]):
                print("--------------------------------------------------------------------------------"):
                print("Checking termination condition for numerator"):
                print("nops(R_num_eval[",k,"]): ", nops(R_num_eval[k])):
                print("R_num_eval[",k,"]: ", R_num_eval[k]);
                print("terms_num_eval[",k,"]: ", terms_num_eval[k]);
                print("T: ", T);
                
                # Add check for empty R_num_eval[k]
                if R_num_eval[k] = [] then
                    print("MRFI: Empty roots list for numerator component ", k):
                    numerator_done[k] := false:
                    next:
                end if:

                # print("R_num_eval[",k,"]: ", R_num_eval[k]);
                # print("R_num_eval[",k,"][1]: ", R_num_eval[k][1]);
                # print("R_num_eval[",k,"][1][1]: ", R_num_eval[k][1][1]);
                if nops(R_num_eval[k]) > 0 and R_num_eval[k][1][1] = 0 then
                    R_num_eval[k] := remove(x->x=[0,1], R_num_eval[k]):
                    terms_num_eval[k] := terms_num_eval[k] - 1:
                end if:
                # if nops(R_num_eval[k]) = terms_num_eval[k] and terms_num_eval[k] < T then

                if nops(R_num_eval[k]) = terms_num_eval[k] and terms_num_eval[k] <= T then
                    print("MRFI Numerator component ",k," recovered!"):
                    numerator_done[k] := true:
                end if:
                # else
                #     all_num_done := false:
            # end if:
            print("____________________________________________________________________________"):
        end do:
        print("===================================================================================="):
        print("MFRI processing denominators now"):
        # Process denominators
        lprint("MRFI den_eval: ",den_eval):
        lambda_den_eval := []:
        terms_den_eval := []:
        R_den_eval := []:
        
        all_den_done := true:
        lprint("MRFI common_den_flag: ",common_den_flag):
        if common_den_flag then
            lprint("MRFI denominator_done: ",denominator_done):
            denominator_done[1] := false:
            if not denominator_done[1] then
                temp := BMEA(den_eval[1], p, Z):
                lprint("MRFI temp (denominator)=",temp):
                lambda_den_eval[k] := [temp]:
                terms_den_eval[k] :=degree(temp, Z):
                R_den_eval[k] := Roots(temp) mod p:
                lprint("MRFI lambda_den_eval: ",lambda_den_eval):
                # terms_den_eval := [degree(temp, Z)]:
                lprint("MRFI terms_den_eval: ",terms_den_eval):
                # R_den_eval := [Roots(temp) mod p]:
                lprint("MRFI R_den_eval: ",R_den_eval):
                
                if nops(R_den_eval[1]) > 0 and R_den_eval[1][1][1] = 0 then
                    R_den_eval[1] := remove(x->x=[0,1], R_den_eval[1]):
                    terms_den_eval[1] := terms_den_eval[1] - 1:
                end if:
                
                if nops(R_den_eval[1]) = terms_den_eval[1] and terms_den_eval[1] < T then
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do:
                # else
                #     all_den_done := false:
                end if:
            end if:
        else
            for k from 1 to num_eqn do
                if  denominator_done[k] then 
                    print("MRFI Skipping BMEA for denominator component ",k," as already done"): 
                    next:
                end if;
                # temp := BMEA(den_eval[k], p, Z):
                # accumulator list
                # lambda_den_eval := [op(lambda_den_eval), temp]:
                # terms_den_eval := [op(terms_den_eval), degree(temp, Z)]:
                # R_den_eval:= [op(R_den_eval), Roots(temp) mod p]:
                lambda_den_eval := [op(lambda_den_eval), BMEA(den_eval[k], p, Z)]:
                terms_den_eval := [op(terms_den_eval), degree(lambda_den_eval[k], Z)]:
                R_den_eval := [op(R_den_eval), Roots(lambda_den_eval[k]) mod p]:

                # hash
                # lambda_den_eval[k] := [temp]:
                # terms_den_eval[k] := degree(temp, Z):
                # R_den_eval[k] := Roots(temp) mod p:

                    print("--------------------------------------------------------------------------------"):
                    print("Checking termination condition for denominator"):
                    print("nops(R_den_eval[",k,"]): ", nops(R_den_eval[k])):
                    print("R_den_eval[",k,"][1][1]: ", R_den_eval[k][1][1]);
                    print("terms_den_eval[",k,"]: ", terms_den_eval[k]);
                    print("T: ", T);

                # Add check for empty R_num_eval[k]
                if R_den_eval[k] = [] then
                    print("MRFI: Empty roots list for numerator component ", k):
                    denominator_done[k] := false:
                    next:
                end if:

                    if nops(R_den_eval[k]) > 0 and R_den_eval[k][1][1] = 0 then
                        R_den_eval[k] := remove(x->x=[0,1], R_den_eval[k]):
                        terms_den_eval[k] := terms_den_eval[k] - 1:
                    end if:
                    
                    if nops(R_den_eval[k]) = terms_den_eval[k] and terms_den_eval[k] < T then
                        print("MRFI Denominator component ",k," recovered!"):
                        denominator_done[k] := true:
                    # else
                    #     all_den_done := false:
                    end if:
                # end if:
            end do:
        end if:

        # for i from 1 to num_eqn do 
        #     if numerator_done[i] = false  then
        #         lambda_num_eval[i]:=[];
        #         terms_num_eval[i]:=[];
        #         R_num_eval[i]:=[];
        #     end if;
        #     if denominator_done[i] = false  then
        #         lambda_den_eval[i]:=[];
        #         terms_den_eval[i]:=[];
        #         R_den_eval[i]:=[];
        #     end if;
        # end do;
        print("--------------------------------------------------------------------------------"):
        for i from 1 to num_eqn do 
            bmea_done[i] := numerator_done[i] and denominator_done[i]:
        end do:
        print("BMEA done status: ",bmea_done):
        for i from 1 to num_eqn do 
            all_done:=all_done and bmea_done[i]:
        end do:
        print("All done status: ",all_done):
        if all_done then
            lprint("MRFI All components recovered!"):
            break:
        end if:
        print("numerator_done: ",numerator_done):
        print("denominator_done: ",denominator_done):
        print("R_num_eval: ",R_num_eval):
        print("R_den_eval: ",R_den_eval):
        print("terms_num_eval: ",terms_num_eval):
        print("terms_den_eval: ",terms_den_eval):
        print("lambda_num_eval: ",lambda_num_eval):
        print("lambda_den_eval: ",lambda_den_eval):
        # break;
        T_old := 2*T:
        T := T*2:
        print("______________________________________________________________________________"):
        if( T > 2^6) then break; end if; # Safety break to avoid infinite loops during testing
    end do:
    
    # Extract roots and build polynomials
    Roots_num_eval := [seq([seq(r[1], r in R_num_eval[k])], k=1..num_eqn)]:
    print("MRFI Roots_num_eval: ",Roots_num_eval):
    
    if common_den_flag then
        Roots_den_eval := [[seq(r[1], r in R_den_eval[1])]]:
        for k from 2 to num_eqn do
            Roots_den_eval := [op(Roots_den_eval), Roots_den_eval[1]]:
        end do:
    else
        Roots_den_eval := [seq([seq(r[1], r in R_den_eval[k])], k=1..num_eqn)]:
    end if:
    print("MRFI Roots_den_eval: ",Roots_den_eval):
    # Generate monomials and coefficients
    num_mono_list := []:
    coeff_num_eval := []:
    final_num_eval := []:
    
    for k from 1 to num_eqn do
        temp := generate_monomials(Roots_num_eval[k], num_vars, Primes, vars):
        if temp = FAIL then return FAIL: end if:
        num_mono_list := [op(num_mono_list), temp]:
        
        coeff_num_eval := [op(coeff_num_eval), 
            Zippel_Vandermonde_solver(num_eval[k], terms_num_eval[k], 
                                     Roots_num_eval[k], lambda_num_eval[k], p)]:
        
        final_num_eval := [op(final_num_eval), 
            construct_final_polynomial(coeff_num_eval[k], num_mono_list[k])]:
    end do:
    
    # Generate denominators
    den_mono_list := []:
    coeff_den_eval := []:
    final_den_eval := []:
    
    if common_den_flag then
        temp := generate_monomials(Roots_den_eval[1], num_vars, Primes, vars):
        if temp = FAIL then return FAIL: end if:
        
        coeff_den := Zippel_Vandermonde_solver(den_eval[1], terms_den_eval[1],
                                               Roots_den_eval[1], lambda_den_eval[1], p):
        
        temp_den := construct_final_polynomial(coeff_den, temp):
        final_den_eval := [seq(temp_den, k=1..num_eqn)]:
    else
        for k from 1 to num_eqn do
            temp := generate_monomials(Roots_den_eval[k], num_vars, Primes, vars):
            if temp = FAIL then return FAIL: end if:
            
            coeff_den := Zippel_Vandermonde_solver(den_eval[k], terms_den_eval[k],
                                                   Roots_den_eval[k], lambda_den_eval[k], p):
            
            final_den_eval := [op(final_den_eval), 
                construct_final_polynomial(coeff_den, temp)]:
        end do:
    end if:
    
    # Normalize
    for k from 1 to num_eqn do
        den_lc := lcoeff(final_den_eval[k], order=grlex(seq(vars[i], i=1..num_vars))):
        den_lc_inv := 1/den_lc mod p:
        final_num_eval[k] := final_num_eval[k] * den_lc_inv mod p:
        final_den_eval[k] := final_den_eval[k] * den_lc_inv mod p:
    end do:
    
    lprint("MRFI ========================================"):
    lprint("MRFI RECOVERY COMPLETE"):
    lprint("MRFI ========================================"):
    
    return final_num_eval, final_den_eval:
end proc:
