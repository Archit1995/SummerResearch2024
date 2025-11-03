
########################################################
# 2. NDSA - For Vector Black Boxes
########################################################

with(LinearAlgebra):

NDSA:=proc(B,sigma_,beta_,num_var,p,num_points)
    local correct_degree,T,alpha,m,phi_,_phi,Y,u,f,g,dq,lcg,np,nv,i,r,
          lin_sys,temp,result,count,M,row,col,DQ:
    lprint("In NDSA");
    correct_degree:=false:
    lin_sys:=false:
    T:=num_points:
    temp:=[]:
    result:=[]:
    count:=0:
    
    while(not(correct_degree)) do
        count:=count+1:
        lprint("T:= ",T):
        r:=rand(p):
        # alpha:=[seq(r()+r() mod p,i=1..T)]:
        alpha:=[seq(i+1 mod p,i=1..T)]:
        lprint("NDSA: alpha: ",alpha):
        m:=expand(product(x-alpha[j],j=1..T)) mod p:
        lprint("NDSA:m: ",m):
        _phi:=projection_image_phi(num_var,alpha,beta_,sigma_,p,T):
        lprint("NDSA:_phi: ",_phi):
        Y:=[seq(B(_phi[i],p),i=1..T)]:
        M:=convert(Y,Matrix):
        # lprint("Matrix M: ",M):
        row,col:=Dimension(M):
        # lprint("row: ",row, " col: ",col):
        
        if row =1 then 
            u:=Interp(alpha,Y,x)mod p:
        else
            lin_sys:=true: 
            u:=get_u(M,col,alpha,p):
        end if:
        lprint("NDSA: u: ",u):
        
        if lin_sys = false then  
            result:=[MQRFR(m,u,0,1,p)]:
            dq:=result[1][3]:
        else 
            for i from 1 to nops(u) do 
                temp:=[op(temp),MQRFR(m,u[i],0,1,p)]:
                result:=[op(result),temp]:
                temp:=[]:
            end do:
            lprint("NDSA:result: ",result):
            DQ:=[seq(result[i][3],i=1..nops(result))]:
            lprint("NDSA:DQ: ",DQ):
            dq:=min(DQ):
            lprint("NDSA:Minimum dq: ",dq):
        end if:
        
        if dq > 1 then  
            lprint("NDSA: Termination condition met"):    
            lprint("NDSA: T:= ",T):
            lprint("NDSA: num_points:= ",num_points):      
            if num_points <> T then  
                
                return result,lin_sys:
            else
                return result,lin_sys: 
            end if:
        else 
            lprint("NDSA:MQRFR failed. Trying again with more points"):
            T:=T*2:
            result:=[]:
            DQ:=[]:
        end if:
        if(count= 4) then break: end if:
    end do:
end proc:


