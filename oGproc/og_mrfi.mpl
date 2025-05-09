MRFI:=proc(B,num_vars::integer,vars::list,p::integer)
    local i,j,Primes,shift_,memo,num,den,num_points_mqrfr,numerator_done,denominator_done,
    T,T_old,j_init,f,g,lg,num_eval,den_eval,deg_num,deg_den,R_num,R_den,r_,ii,lambda_num,lambda_den
    ,terms_num,terms_den,Roots_num,Roots_den,num_mono,den_mono,coeff_num,coeff_den,den_lc,r,den_lc_inv:
    print("In MRFI"):
    r_:=rand(p):
    Primes:=generate_primes(1,num_vars):
    shift_:=[seq(r_(),i=1..num_vars-1)]:
    memo:=[]:
    num:=[]:
    den:=[]:
    num_points_mqrfr:=0:
    numerator_done:=false:
    denominator_done:=false:
    T:=4:
    T_old:=0:
    j_init:=2:
    f,g,lg,num_points_mqrfr:=NDSA(B,[seq(1,i=1..num_var)],shift_,num_var,p,T): 
    num_eval:=eval(f,x=1):
    num:=[op(num),num_eval]:
    den_eval:=eval(g,x=1):
    den:=[op(den),den_eval]:
    deg_num:=degree(f,x):
    deg_den:=degree(g,x):
    num_points_mqrfr:=deg_num+deg_den+2:
    print("num_points_mqrfr = ",num_points_mqrfr):
    print("den= ",den):
    print("num= ",num):
    while true do
        print("T=",T):
        print("T_old=",T_old):  
        # print("Primess: ",Primes):
        # [2^j,3^j]
        # for ii from t_old to 2*t do 
        for ii from T_old to 2*T-1 do 
        memo:=[op(memo),[seq(Primes[j]^ii,j=1..num_var)]] mod p:
        end do;
        print("memo=",memo):
        # break;
        # Primess:=[op(Primess),[seq(modp(Primes[ii]^i,p),ii=1..num_var)]]:
        print("numelems(memo) = ",numelems(memo));
        # get degree here 
        
        print("j_init=",j_init):
    
        for j from j_init to numelems(memo) do 
            print("Back in get_lambda j=",j):
            # f,g,lg:=evaluation_num_den(B,memo[j][1],shift_,num_var,p):
            print("memo[",j,"]=",memo[j]):
            f,g,lg:=NDSA(B,memo[j],shift_,num_var,p,num_points_mqrfr): 
            if not(numerator_done) then 
                # print("numerator_done= ",numerator_done):
                num_eval:=eval(f,x=memo[j][1]):
                # print("num_eval=",num_eval):
                num:=[op(num),num_eval]:
                print("num= ",num):
            end if:
            
            if not(denominator_done) then 
                # print("denominator_done= ",denominator_done):
                den_eval:=eval(g,x=memo[j][1]):
                # print("den_eval=",den_eval):
                den:=[op(den),den_eval]:               
                print("den= ",den):
            end if:

            print("________________________________________________________________"):
        end do:

     
        # Call BMEA only on the polynomial which caused the done, while keeping the value of the other polynomial.
        if not(numerator_done) then 
            lambda_num := BMEA(num,p,Z):
            terms_num:=degree(lambda_num,Z):
            # terms_num[i]:=degree(lambda_num,Z):
            R_num := Roots(lambda_num) mod p:
            print("lambda_num= ",lambda_num):
            print("terms_num[i]= ",terms_num):
            print("R_num= ",R_num):
        end if:
        if not(denominator_done) then 
            lambda_den := BMEA(den,p,Z):
            print("lambda_den= ",lambda_den): 
            terms_den:=degree(lambda_den,Z):
            # terms_den[i]:=degree(lambda_den,Z):
            print("terms_den[i]= ",terms_den):
            # R_num := Roots(lambda_num) mod p:
            R_den := Roots(lambda_den) mod p:
            print("R_den= ",R_den):
        end if:

        print("num of roots of lambda_num =",nops(R_num)):
        print("num of roots of lambda_den =",nops(R_den)):
        if nops(R_num)=terms_num and terms_num < T then 
            print("Numerator Success!!"):
            numerator_done:=true:
        end if:
        # if nops(R_den)=terms_den[i] then 
        if nops(R_den)=terms_den and terms_den < T  then 
            print("Denominator Success!!"):
            denominator_done:=true:
        end if:
        if numerator_done and denominator_done then 
            print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");  
            # print("terms_num= ",terms_num):
            # print("terms_den= ",terms_den):
            # print("lambda_num= ",lambda_num):
            # print("lambda_den= ",lambda_den):
            # print("R_num= ",R_num):
            # print("R_den= ",R_den):
            # print("num= ",num):
            # print("den= ",den):
            # return terms_num,terms_den,lambda_num,lambda_den,R_num,R_den,num,den:
            break:
        end if:
        T_old:=2*T:
        j_init:=T_old+1:
        T:=T*2:   
        print("T=",T):
        print("T_old=",T_old):
    end do:
    Roots_num := [ seq(r[1], r in R_num ) ]:
    Roots_den := [ seq(r[1], r in R_den ) ]:
    
    num_mono:= generate_monomials(Roots_num,num_var,Primes,vars):
    print("num_mono: ",num_mono):
    den_mono:= generate_monomials(Roots_den,num_var,Primes,vars):
    print("den_mono: ",den_mono):
    if num_mono = FAIL or den_mono = FAIL then 
        return FAIL:
    end if:
    print("terms_num: ",terms_num):
    # print(type(Roots_num,list)):
    print("Roots_num: ",Roots_num):
    coeff_num:= Zippel_Vandermonde_solver(num,terms_num,Roots_num,lambda_num,p):
    print("coeff_num: ",coeff_num):
    coeff_den:= Zippel_Vandermonde_solver(den,terms_den,Roots_den,lambda_den,p):
    print("coeff_den: ",coeff_den):
    num:= construct_final_polynomial(coeff_num,num_mono):
    den:= construct_final_polynomial(coeff_den,den_mono):
    # print("num: ",num):
    # print("den: ",den):
    den_lc:=lcoeff(den,order=grlex(x,y,z)):
    den_lc_inv:=1/den_lc mod p:
    # print("den_lc: ",den_lc):
    # print("num= ",num*den_lc_inv mod p):
    # print("den= ",den*den_lc_inv mod p):
    return num*den_lc_inv mod p, den*den_lc_inv mod p:
end proc:
