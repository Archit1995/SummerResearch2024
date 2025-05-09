# 2. Generating a prime for each variable 
generate_evaulation_primes:=proc(init_,n)
    print("In generate_evaulation_primes"):
    local sigma_,m,i:
    m:=init_:

    # print("m=",m):
    sigma_:=Vector(n,0):
    for i from 1 to n do 
        sigma_[i]:=nextprime(m):
        # print("sigma_[",i,"]=",sigma_[i]):
        m:=sigma_[i]:
    end do:
    print("sigma_=",sigma_):    
return convert(sigma_,list),m:

end proc: