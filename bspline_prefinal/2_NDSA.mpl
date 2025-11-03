
########################################################
# 2. NDSA - For Vector Black Boxes
########################################################

with(LinearAlgebra):

NDSA:=proc(B,sigma_,beta_,num_var,p,num_points,num_eqn)
    local correct_degree,T,alpha,m,phi_,_phi,Y,u,f,g,dq,lcg,np,nv,i,r,
          lin_sys,temp,result,count,M,row,col,DQ,MQRFR_done:
    print("In NDSA");
    # print("NDSA: num_eqn:= ",num_eqn);
    MQRFR_done:=[seq(false, i=1..num_eqn)]:
    correct_degree:=false:
    lin_sys:=false:
    T:=num_points:
    temp:=[]:
    # result:=[]:
    result:=[seq([], i=1..num_eqn)]:
    # result:=Array(1..num_eqn,[]):
    count:=0:
    
    while(not(correct_degree)) do
        count:=count+1:
        
        lprint("T:= ",T):
        r:=rand(p):
        # r:=rand(2^31-1):
        alpha:=[seq(r() mod p,i=1..T)]:
        # alpha:=[seq(i+1 mod p,i=1..T)]:
        # lprint("NDSA: alpha: ",alpha):
        m:=expand(product(x-alpha[j],j=1..T)) mod p:
        # lprint("NDSA:m: ",m):
        _phi:=projection_image_phi(num_var,alpha,beta_,sigma_,p,T):
        # lprint("NDSA:_phi: ",_phi):
        Y:=[seq(B(_phi[i],p),i=1..T)]:
        # print("NDSA: Y = ",Y);
        # M:=convert(Y,Matrix):
        M:=Matrix(Y):
        # lprint("Matrix M: ",M):
        row,col:=Dimension(M):
        # lprint("row: ",row, " col: ",col):
        
        if row =1 then 
            u:=Interp(alpha,Y,x)mod p:
        else
            lin_sys:=true: 
            u:=get_u(M,col,alpha,p):
        end if:
        # lprint("NDSA: u: ",u):
        # print("nops u: ",nops(u));  
        
        if lin_sys = false then  
            result:=[MQRFR(m,u,0,1,p)]:
            dq:=result[1][3]:
        else 
            for i from 1 to nops(u) do 
            # put a check for MQRFR_done
                if(MQRFR_done[i]) then 
                    lprint("NDSA: Skipping MQRFR for equation ",i," as already done"):
                    next:
                end if;
                temp:=[op(temp),MQRFR(m,u[i],0,1,p)]:
                # print("NDSA: result temp before assignment: ",result):
                # print("NDSA: temp",i,": ",temp):
                # result:=[op(result),temp]:
                result[i]:=temp:
                # print("NDSA: result",i,": ",result[i]):
                temp:=[]:
            end do:
            # i:=1:
            # print("i=",i):
            # lprint("NDSA:entire result: ",result):
            # result_list := convert(result,list):

            DQ:=[seq(result[i][3],i=1..nops(result))]:
            for i from 1 to numelems(DQ) do 
                if DQ[i] > 1 then 
                    MQRFR_done[i]:=true:
                    print("NDSA: MQRFR successful for equation ",i):
                end if;
            end do:
            # lprint("NDSA:DQ: ",DQ):
            dq:=min(DQ):
            # lprint("NDSA:Minimum dq: ",dq):
        end if:
        
        if dq > 1 then  
            print("NDSA: Termination condition met"):       
            return result,lin_sys:
            # if num_points <> T then  
            #     return result,T,lin_sys:
            # else
            #     return result,lin_sys: 
            # end if:
        else 
            print("NDSA:MQRFR failed. Trying again with more points"):
            T:=T*2:
            # result:=[]:
            print("NDSA: mqrfr_done before reset: ",MQRFR_done):
            for i from 1 to num_eqn do 
                if MQRFR_done[i]= false then 
                    lprint("NDSA: Resetting result for equation ",i):
                    result[i]:=[]:
                end if;
            end do;
            print("_______________________________________________"):
            # result:=[seq([], i=1..num_eqn)]:
            DQ:=[]:
        end if:
        # if(count= 4) then break: end if:
    end do:
end proc:


