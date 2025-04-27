# read "./projection_phi.mpl":
# anchor points=[[2,3],[4,9],[8,27],[16,81],[32,243]...]:
# Shift=[5]
evaluation_num_den:=proc(B,num_points,p)
    local correct_degree,T,alpha,m,phi_,_phi,Y,u,f,g,dq,lcg,np,nv:
    print("---------------------------------------------------------------------");
    print("In evaluation_num_den");

    correct_degree:=false:
    T:=num_points:
    while(not(correct_degree)) do
        
        print("T:= ",T):
# we need to come up with some values of alpha
        # alpha:=[seq(i-1,i=1..T)]:
        r:=rand(p):
        alpha:=[seq(r()mod p,i=1..T)]:
        print("alpha: ",alpha):
        print("alpha[1]: ",alpha[1]):
        m:=expand(product(x-alpha[i],i=1..T)) mod p:
        print("m: ",m):

        # _phi:=projection_image_phi(num_var,alpha,shift_,anchor_points,p,T):
        print("_phi: ",_phi):
        Y:=[seq(B(alpha[i],p),i=1..T)]:
        print("Y: ",Y):
        u:=Interp(alpha,Y,x)mod p:
        print("u: ",u):
        # u:=u/lcoeff(u) mod p:
        # m:=m*(1/lcoeff(u)) mod p:
        f,g,dq,lcg:=MQRFR(m,u,0,1,p):
        # f,g,dq,lcg:=MQRFR(m,u,0,1/lcoeff(u) mod p,p):
        print("dq: ",dq):
        if dq > 1 then 
            break:
        else 
            print("MQRFR failed. Trying again with more points"):
            T:=T*2:
        end if:
        print("====================================================="):
    end do:
    # if num_points <> T then  
    #     return f,g,lcg,T:
    #     # return m,u,f,g,lcg,T:

    # else
    return f,g,lcg:
    # Checking for Ratrecon
        # return m,u,f,g,lcg:
    # end if:
end proc:
    # We need to calculate m, y and u  