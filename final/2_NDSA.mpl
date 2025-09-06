NDSA:=proc(B,sigma_,beta_,num_var,p,num_points)
    local correct_degree,T,alpha,m,phi_,_phi,Y,u,f,g,dq,lcg,np,nv,i,r:
    # print("---------------------------------------------------------------------");
    # print("In evaluation_num_den");
    correct_degree:=false:
    T:=num_points:

    while(not(correct_degree)) do
        # print("T:= ",T):
        r:=rand(p):
        alpha:=[seq(r()+r() mod p,i=1..T)]:
        # print("alpha: ",alpha):
        m:=expand(product(x-alpha[i],i=1..T)) mod p:
        # print("m: ",m):
        _phi:=projection_image_phi(num_var,alpha,beta_,sigma_,p,T):
        print("_phi: ",_phi):
        Y:=[seq(B(_phi[i],p),i=1..T)]:
        u:=Interp(alpha,Y,x)mod p:
        f,g,dq,lcg:=MQRFR(m,u,0,1,p):
        print("f: ",f):
        print("g: ",g):
        print("dq: ",dq):
        if dq > 1 then            
            correct_degree:=true:
            break:
        else 
            print("MQRFR failed. Trying again with more points"):
            T:=T*2:
        end if:
        # print("====================================================="):
    end do:
    if num_points <> T or correct_degree = false then  
    # print("====================================================="):
        return f,g,lcg,T:
    else
    # print("====================================================="):
    return f,g,lcg:
    end if:
end proc:
