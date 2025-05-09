
get_sigma:=proc(B,num_var,p)
    local sigma_init,eval_sigma,sigma_,i,init_:
    init_:=1:
    denominator_zero:=true:
    print("In get_sigma");
    sigma_,next_:=generate_evaulation_primes(init_,num_var):
    init_:=next_:
    print("init_=",init_):
    print("sigma_=",sigma_):
    while denominator_zero do 
    try 
        denominator_zero:=false:
        eval_sigma:=B(sigma_,p):
    catch "Denominator is zero":
        denominator_zero:=true:
        sigma_,next_:=generate_evaulation_primes(init_,num_var):
        init_:=next_:
        print("Denominator is zero"):
        print("sigma_=",sigma_):
        # return 0:
    # catch:
    #     print("default case"):
        
    end try:
    #
    end do:
    
    return sigma_:
end proc: