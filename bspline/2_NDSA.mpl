with(LinearAlgebra):
NDSA:=proc(B,sigma_,beta_,num_var,p,num_points)
    local correct_degree,T,alpha,m,phi_,_phi,Y,u,f,g,dq,lcg,np,nv,i,r:
    # print("---------------------------------------------------------------------");
    print("In NDSA");
    correct_degree:=false:
    lin_sys:=false:
    T:=num_points:
    # result and DQ are buffers. 
    temp:=[]:
    result:=[]:
    count:=0:
    while(not(correct_degree)) do
        count:=count+1:
        print("num_points:= ",num_points):
        print("T:= ",T):
        r:=rand(p):
        alpha:=[seq(r()+r() mod p,i=1..T)]:
        print("alpha: ",alpha):
        m:=expand(product(x-alpha[j],j=1..T)) mod p:
        print("m: ",m):
        _phi:=projection_image_phi(num_var,alpha,beta_,sigma_,p,T):
        print("_phi: ",_phi):
        Y:=[seq(B(_phi[i],p),i=1..T)]:
        M:=convert(Y,Matrix):
       
        row,col:=Dimension(M):
        print("Y: ",Y):
        if row =1 then 
            u:=Interp(alpha,Y,x)mod p:
        else
            lin_sys:=true: 
            u:=get_u(M,col,alpha,p):
        end if:
        print("u: ",u):
        # print(eval(B)):
        # print(op(4,eval(B))):
        
        if lin_sys = false then  
            # f,g,dq,lcg:=MQRFR(m,u,0,1,p):
            result:=[MQRFR(m,u,0,1,p)]:
        else 
            for i from 1 to nops(u) do 
                temp:=[op(temp),MQRFR(m,u[i],0,1,p)]:
                result:=[op(result),temp]:
                temp:=[]:
            end do:
            print("result=",result):
            DQ:=[seq(result[i,3],i=1..nops(result))]:
            print("DQ: ",DQ):
            dq:=min(DQ):
            print("dq: ",dq):

        end if:
        # print("result: ",result):
        if  dq > 1 then  
            print("Termination condition met"):          
            if num_points <> T then  
                print("in num_points <> T"):
                # print("Returning: result=", result, " T=", T, " lin_sys=", lin_sys):
                return result,T,lin_sys:
            else
                return result: 
            end if:
        else 
            print("MQRFR failed. Trying again with more points"):
            T:=T*2:
            # Resetting buffers result and DQ
            result:=[]:
            DQ:=[]:
        end if:
        print("====================================================="):
        if(count= 4) then break: end if:
    end do:
     
end proc:
