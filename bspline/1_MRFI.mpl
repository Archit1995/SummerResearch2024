MRFI:=proc(B,num_vars::integer,num_eqn::integer,vars::list,p::integer)
    local i,j,Primes,shift_,sigma_,num,den,num_points_mqrfr,numerator_done,denominator_done,
    T,T_old,j_init,f,g,lg,num_eval,den_eval,deg_num,deg_den,R_num,R_den,r_,ii
    ,terms_num,terms_den,Roots_num,Roots_den,num_mono,den_mono,coeff_num,coeff_den,den_lc,r,den_lc_inv:
    print("In MRFI"):
    lambda_num:=[]:
    lambda_den:=[]:
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
    # numerator_done:=[seq(false,i=1..num_eqn)]:
    # print("numerator_done: ",numerator_done):
    denominator_done:=false:
    # denominator_done:=[seq(false,i=1..num_vars)]:
    
    T:=4:
    # T_old:=0:
    T_old:=1:
    # print("Yuhuuu", NDSA(B,[seq(1,i=1..num_vars)],shift_,num_vars,p,T)):
    # mqrfr_results=[num_poly(x),denom_poly(x),qmax,leading denom_coeff] is a buffer: Shoudl it not be qmin?

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
        
        # F:=[seq(mqrfr_results[i][1],i=1..nops(mqrfr_results))]:

        temp:=[]:
        temp:=[seq(mqrfr_results[i][1],i=1..nops(mqrfr_results))]:
        F:=[op(F),temp]:
        temp:=[]:
        # G:=[seq(mqrfr_results[i][2],i=1..nops(mqrfr_results))]:
        temp:=[seq(mqrfr_results[i][2],i=1..nops(mqrfr_results))]:
        G:=[op(G),temp]:
        temp:=[]:
    
        # lg:=[seq(mqrfr_results[i][3],i=1..nops(mqrfr_results))]:
        temp:=[seq(mqrfr_results[i][3],i=1..nops(mqrfr_results))]:
        lg:=[op(lg),temp]:
        temp:=[]:
        print("F: ",F):
        print("G: ",G):
        # Collect the degrees of the numerators and denominators
        deg_num:=[seq(degree(F[1][i],x),i=1..nops(F))]:

        deg_den:=[seq(degree(G[1][i],x),i=1...nops(G))]:
        print("deg_num: ",deg_num):
        print("deg_den: ",deg_den):
        # Calculate the number of points needed for MQRFR
        # num_points_mqrfr:=[seq(deg_num[i]+deg_den[i]+2,i=1..nops(deg_num))]:
        num_points_mqrfr:=max(deg_num)+max(deg_den)+2:
        print("num_points_mqrfr: ",num_points_mqrfr):
        # Evaluate the numerators and denominators at x=1
        num_eval:=[seq(eval(F[i],x=1),i=1..nops(F))]:
        
        den_eval:=[seq(eval(G[i],x=1),i=1..nops(G))]:
        print("num_eval: ",num_eval):
        print("den_eval: ",den_eval):
        # F,G,lg,deg_num,deg_den,num_points_mqrfr,num_eval,den_eval:=update_buffers(mqrfr_results):
        #refreshing buffer 
        mqrfr_results:=[]:
    end if:
   count:=0:
    while true do
        print("T=",T):
        print("T_old=",T_old):  

        for j from T_old to 2*T-1 do   #powers of prime [2^j,3^j,...,p^j]

            sigma_:=[op(sigma_),[seq(Primes[i]^j mod p,i=1..nops(Primes))]]:
            print("sigma_[",j,"]=",sigma_[j]):
            print("F: ",F):
            print("nops(F)=",nops(F)):
        #     # f,g,lg:=NDSA(B,sigma_[j],shift_,num_vars,p,num_points_mqrfr): 
            # for k from 1 to nops(num_points_mqrfr) do 
                # print("k=",k):
                mqrfr_results:=NDSA(B,sigma_[j],shift_,num_vars,p,num_points_mqrfr): 
                print("mqrfr_results: ",mqrfr_results):
                temp:=[]:
                temp:=[seq(mqrfr_results[i][1],i=1..nops(mqrfr_results))]:
                F:=[op(F),temp]:
                # F:=[temp,nops(F)]:
                temp:=[]:
                print("F: ",F):
            # F:=[seq(mqrfr_results[i][1],i=1..nops(mqrfr_results))]:
                temp:=[seq(mqrfr_results[i][2],i=1..nops(mqrfr_results))]:
                G:=[op(G),temp]:
                temp:=[]:

            # G:=[seq(mqrfr_results[i][2],i=1..nops(mqrfr_results))]:
                print("G: ",G):
                temp:=[seq(mqrfr_results[i][3],i=1..nops(mqrfr_results))]:
                lg:=[op(lg),temp]:
                temp:=[]:
            # lg:=[seq(mqrfr_results[i][3],i=1..nops(mqrfr_results))]:
            # print("F: ",F):
            # print("G: ",G):
            # # Collect the degrees of the numerators and denominators
            # deg_num:=[seq(degree(F[i],x),i=1..nops(F))]:
            # deg_den:=[seq(degree(G[i],x),i=1...nops(G))]:
            # print("deg_num: ",deg_num):
            # print("deg_den: ",deg_den):
            # # Calculate the number of points needed for MQRFR
            # # num_points_mqrfr:=[seq(deg_num[i]+deg_den[i]+2,i=1..nops(deg_num))]:
            # num_points_mqrfr:=max(deg_num)+max(deg_den)+2:
            # print("num_points_mqrfr: ",num_points_mqrfr):
            # # Evaluate the numerators and denominators at x=1
                temp:=[seq(eval(F[j+1][i],x=1),i=1..nops(F[j+1]))]:
                num_eval:=[op(num_eval),temp]:
                temp:=[]:
                print("num_eval: ",num_eval):
                temp:=[seq(eval(G[j+1][i],x=1),i=1..nops(G[j+1]))]:
                den_eval:=[op(den_eval),temp]:
                temp:=[]:
                print("den_eval: ",den_eval):

                num:=[seq(num_eval[..,i],i=1..nops(num_eval[1]))]:
                den:=[seq(den_eval[..,i],i=1..nops(den_eval[1]))]:
                print("num: ",num):
                print("den: ",den):
            # num:=[op(num),num_eval]:;
                # num_eval:=[[seq(eval(F[i],x=1),i=1..nops(F))]]:
            # den_eval:=[[seq(eval(G[i],x=1),i=1..nops(G))]]:
            # print("num_eval: ",num_eval):
            # print("den_eval: ",den_eval):
            # mqrfr_results:=[]:
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
        if not(numerator_done) then 
            for i from 1 to nops(num) do 
                lambda_num:=[op(lambda_num), BMEA(num[i],p,Z)]:
            end do;
            print("lambda_num: ",lambda_num):
            terms_num:=[seq(degree(lambda_num[i],Z),i=1..nops(lambda_num))]:
            print("terms_num: ",terms_num):
            R_num:=[seq(Roots(lambda_num[i]) mod p,i=1..nops(lambda_num))]:
            print("R_num: ",R_num):
            num_roots_num:=[seq(nops(R_num[i]),i=1..nops(R_num))]:
            print("num_roots_num: ",num_roots_num):
             print("terms_num: ",nops(terms_num)):
            for ii from 1 to nops(terms_num) do 
                if num_roots_num[ii]=terms_num[ii] and terms_num[ii] < T then 
                    # print("Numerator Success!! for i=",ii):
                    numerator_done:=true:
                else 
                    
                    numerator_done:=false:
                    print("BEMA failed for numerator for i="):
                    lambda_num:=[]:
                    terms_num:=[]:
                    R_num:=[]:
                    num_roots_num:=[]:
                    break:
                
                end if;
            end do;
            print("numerator_done: ",numerator_done):
        end if;
        if not(denominator_done) then
            for i from 1 to nops(den) do 
                lambda_den:=[op(lambda_den), BMEA(den[i],p,Z)]:
            end do;
            print("lambda_den: ",lambda_den):
            terms_den:=[seq(degree(lambda_den[i],Z),i=1..nops(lambda_den))]:
            print("terms_den: ",terms_den):
            R_den:=[seq(Roots(lambda_den[i]) mod p,i=1..nops(lambda_den))]:
            print("R_den: ",R_den):
            num_roots_den:=[seq(nops(R_den[i]),i=1..nops(R_den))]:
            print("num_roots_den: ",num_roots_den):
            print("terms_den: ",nops(terms_den)):
            for ii from 1 to nops(terms_den) do 
                if num_roots_den[ii]=terms_den[ii] and terms_den[ii] < T then 
                    # print("Denominator Success!! for i=",ii):
                    denominator_done:=true:
                else 
                    print("BEMA failed for denominator for i="):
                    denominator_done:=false:
                    lambda_den:=[]:
                    terms_den:=[]:
                    R_den:=[]:
                    num_roots_den:=[]:
                    break:
                end if;
            end do;
            print("denominator_done: ",denominator_done):
        end if;
        if numerator_done and denominator_done then 
            print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");
            break:
        end if:

        
            
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
        T_old:=2*T:
        T:=T*2:  
        count:=count+1:
    
    # if  numerator_done then break end if; # To avoid infinite loop during testing.
    end do:
    # Roots_num := [ seq(r[1], r in R_num ) ]:
    # Roots_num:=[seq([seq(r[1],r in seq(R_num[i],i=1..nops(R_num))[j])],j=1..nops(R_num))];
    Roots_num:=[[ 2, 3], [ 2, 3, 4, 6, 9], [ 2, 3]];\\

    

    print("Roots_num: ",Roots_num):
    print("vars: ",vars):
    # Roots_den := [ seq(r[1], r in R_den ) ]:

        num_mono:=[];
    for i from 1 to nops(Roots_num)do 
        print("i=",i):
        num_mono:=[op(num_mono),generate_monomials(Roots_num[i],num_var,Primes,vars)]:
    end do:
    print("num_mono: ",num_mono):
    den_mono:=[];
    for i from 1 to nops(Roots_den)do 
        print("i=",i):
        den_mono:=[op(den_mono),generate_monomials(Roots_den[i],num_var,Primes,vars)]:
    end do:
    print("den_mono: ",den_mono):
    if num_mono = FAIL or den_mono = FAIL then 
        return FAIL:
    end if:
    # den_mono:= generate_monomials(Roots_den,num_var,Primes,vars):
    # print("den_mono: ",den_mono):
    # if num_mono = FAIL or den_mono = FAIL then 
    #     return FAIL:
    # end if:
    # print("terms_num: ",terms_num):
    # print("Roots_num: ",Roots_num):
    coeff_num:=[];
    for i from 1 to nops(num)do
        coeff_num:= [op(coeff_num),Zippel_Vandermonde_solver(num[i],terms_num[i],Roots_num[i],lambda_num[i],p)]:
    end do:
    print("coeff_num: ",coeff_num):
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
