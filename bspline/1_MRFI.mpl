MRFI:=proc(B,num_vars::integer,vars::list,p::integer)
    local i,j,Primes,shift_,sigma_,num,den,num_points_mqrfr,numerator_done,denominator_done,
    T,T_old,j_init,f,g,lg,num_eval,den_eval,deg_num,deg_den,R_num,R_den,r_,ii,lambda_num,lambda_den
    ,terms_num,terms_den,Roots_num,Roots_den,num_mono,den_mono,coeff_num,coeff_den,den_lc,r,den_lc_inv:
    print("In MRFI"):
    r_:=rand(p):
    Primes:=[seq(ithprime(i),i=1..num_vars)]:
    shift_:=[seq(r_(),i=1..num_vars-1)]:
    sigma_:=[]:
    num:=[]:
    den:=[]:
    F:=[]:
    G:=[]:
    num_points_mqrfr:=0:
    numerator_done:=false:
    denominator_done:=false:
    
    T:=4:
    # T_old:=0:
    T_old:=1:
    # print("Yuhuuu", NDSA(B,[seq(1,i=1..num_vars)],shift_,num_vars,p,T)):


    mqrfr_results,num_points_mqrfr,lin_sys:=NDSA(B,[seq(1,i=1..num_vars)],shift_,num_vars,p,T):
    print("lin_sys: ",lin_sys):
    print("nops mqrfr_results: ",nops(mqrfr_results)):


    if lin_sys = false then 
        f:=mqrfr_results[1]:
        g:=mqrfr_results[2]:
        lg:=mqrfr_results[3]:
        deg_num:=degree(f,x):
        deg_den:=degree(g,x):
        num_points_mqrfr:=deg_num+deg_den+2:
        num_eval:=eval(f,x=1):
        num:=[op(num),num_eval]:
        den_eval:=eval(g,x=1):
        den:=[op(den),den_eval]:
    else
        print("in lin_sys = true"): 
        # Collect numerators and denominators from the results of MQRFR
        F:=[seq(mqrfr_results[i][1],i=1..nops(mqrfr_results))]:
        G:=[seq(mqrfr_results[i][2],i=1..nops(mqrfr_results))]:
        lg:=[seq(mqrfr_results[i][3],i=1..nops(mqrfr_results))]:
        print("F: ",F):
        print("G: ",G):
        # Collect the degrees of the numerators and denominators
        deg_num:=[seq(degree(F[i],x),i=1..nops(F))]:
        deg_den:=[seq(degree(G[i],x),i=1...nops(G))]:
        print("deg_num: ",deg_num):
        print("deg_den: ",deg_den):
        # Calculate the number of points needed for MQRFR
        # num_points_mqrfr:=[seq(deg_num[i]+deg_den[i]+2,i=1..nops(deg_num))]:
        num_points_mqrfr:=max(deg_num)+max(deg_den)+2:
        print("num_points_mqrfr: ",num_points_mqrfr):
        # Evaluate the numerators and denominators at x=1
        num_eval:=[[seq(eval(F[i],x=1),i=1..nops(F))]]:
        den_eval:=[[seq(eval(G[i],x=1),i=1..nops(G))]]:
        print("num_eval: ",num_eval):
        print("den_eval: ",den_eval):
        #refreshing buffer 
        mqrfr_results:=[]:
    end if:
   
    # while true do
        print("T=",T):
        print("T_old=",T_old):  
        for j from T_old to 2*T-1 do   #powers of prime [2^j,3^j,...,p^j]
            sigma_:=[op(sigma_),[seq(Primes[i]^j mod p,i=1..nops(Primes))]]:
            print("sigma_[",j,"]=",sigma_[j]):
        #     # f,g,lg:=NDSA(B,sigma_[j],shift_,num_vars,p,num_points_mqrfr): 
            # for k from 1 to nops(num_points_mqrfr) do 
                # print("k=",k):
                mqrfr_results:=NDSA(B,sigma_[j],shift_,num_vars,p,num_points_mqrfr): 
                print("mqrfr_results: ",mqrfr_results):
                
            # end do:

            # f:=mqrfr_results[1]:
            # g:=mqrfr_results[2]:
            # lg:=mqrfr_results[3]:
            # if not(numerator_done) then 
            #     num_eval:=eval(f,x=sigma_[j][1]):
            #     num:=[op(num),num_eval]:
            #     # print("num= ",num):
            # end if:
            # if not(denominator_done) then 
            #     den_eval:=eval(g,x=sigma_[j][1]):
            #     den:=[op(den),den_eval]:               
            #     # print("den= ",den):
            # end if:

            # print("________________________________________________________________"):
        end do:
    #         sigma_:=[op(sigma_),[seq(Primes[i]^j mod p,i=1..nops(Primes))]]:
    #         print("sigma_[",j,"]=",sigma_[j]):
    # #         # f,g,lg:=NDSA(B,sigma_[j],shift_,num_vars,p,num_points_mqrfr): 
    #         mqrfr_results,lin_sys:=NDSA(B,sigma_[j],shift_,num_vars,p,num_points_mqrfr): 
    #         print("mqrfr_results: ",mqrfr_results):
    #         print("lin_sys: ",lin_sys):

    #         f:=mqrfr_results[1]:
    #         g:=mqrfr_results[2]:
    #         lg:=mqrfr_results[3]:
    #         if not(numerator_done) then 
    #             num_eval:=eval(f,x=sigma_[j][1]):
    #             num:=[op(num),num_eval]:
    #             # print("num= ",num):
    #         end if:
    #         if not(denominator_done) then 
    #             den_eval:=eval(g,x=sigma_[j][1]):
    #             den:=[op(den),den_eval]:               
    #             # print("den= ",den):
    #         end if:

    #         print("________________________________________________________________"):
        # end do:

    #     # Call BMEA only on the polynomial which caused the done, while keeping the value of the other polynomial.
    #     if not(numerator_done) then 
    #         lambda_num := BMEA(num,p,Z):
    #         terms_num:=degree(lambda_num,Z):
    #         R_num := Roots(lambda_num) mod p:
    #         # print("lambda_num= ",lambda_num):
    #         print("terms_num[i]= ",terms_num):
    #         print("R_num= ",R_num):
    #     end if:
    #     if not(denominator_done) then 
    #         lambda_den := BMEA(den,p,Z):
    #         # print("lambda_den= ",lambda_den): 
    #         terms_den:=degree(lambda_den,Z):
    #         print("terms_den[i]= ",terms_den):
    #         R_den := Roots(lambda_den) mod p:
    #         print("R_den= ",R_den):
    #     end if:

    #     print("num of roots of lambda_num =",nops(R_num)):
    #     print("num of roots of lambda_den =",nops(R_den)):
    #     if nops(R_num)=terms_num and terms_num < T then 
    #         print("Numerator Success!!"):
    #         numerator_done:=true:
    #     end if:
    #     if nops(R_den)=terms_den and terms_den < T  then 
    #         print("Denominator Success!!"):
    #         denominator_done:=true:
    #     end if:
    #     if numerator_done and denominator_done then 
    #         print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");  

    #         break:
    #     end if: 
    #     T_old:=2*T:
    #     T:=T*2:  
    
    # end do:
    # Roots_num := [ seq(r[1], r in R_num ) ]:
    # Roots_den := [ seq(r[1], r in R_den ) ]:
    
    # num_mono:= generate_monomials(Roots_num,num_var,Primes,vars):
    # print("num_mono: ",num_mono):
    # den_mono:= generate_monomials(Roots_den,num_var,Primes,vars):
    # print("den_mono: ",den_mono):
    # if num_mono = FAIL or den_mono = FAIL then 
    #     return FAIL:
    # end if:
    # print("terms_num: ",terms_num):
    # print("Roots_num: ",Roots_num):
    # coeff_num:= Zippel_Vandermonde_solver(num,terms_num,Roots_num,lambda_num,p):
    # print("coeff_num: ",coeff_num):
    # coeff_den:= Zippel_Vandermonde_solver(den,terms_den,Roots_den,lambda_den,p):
    # print("coeff_den: ",coeff_den):
    # num:= construct_final_polynomial(coeff_num,num_mono):
    # den:= construct_final_polynomial(coeff_den,den_mono):
    # den_lc:=lcoeff(den,order=grlex(seq(vars[i],i=1..num_vars))):
    # den_lc_inv:=1/den_lc mod p:
    # print("den_lc= ",den_lc):
    # print("den_lc_inv= ",den_lc_inv):
    # # return num*den_lc_inv mod p, den*den_lc_inv mod p:
    # return num,den 
end proc:
