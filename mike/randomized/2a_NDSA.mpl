########################################################
# 2. NDSA - For Early stopping
########################################################

with(LinearAlgebra):

NDSA:=proc(B,sigma_,beta_,num_var,p,num_points,num_eqn)
    local correct_degree,T,alpha,m,phi_,Psi_alpha,Y,u,f,g,dq,lcg,np,nv,i,r,
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
    count:=0:
    while(not(correct_degree)) do
        count:=count+1:
        
        print("T:= ",T):
        r:=rand(p):
        # r:=rand(2^31-1):
        alpha:=[seq(r() mod p,i=1..T)]:
        # alpha:=[seq(i mod p,i=1..T)]:# for example in Ch4: test case 1
        # alpha:=[seq(i ,i=2..T+1)]:
        lprint("NDSA: alpha: ",alpha):
        m:=expand(product(x-alpha[j],j=1..T)) mod p:
        lprint("NDSA:m: ",m):
        Psi_alpha:=get_point_on_affine_line(num_var,alpha,beta_,sigma_,p,T):
        lprint("NDSA:Psi_alpha: ",Psi_alpha):
        Y:=[seq(B(Psi_alpha[i],p),i=1..T)]:
        print("NDSA: Y = ",Y);
        # M:=convert(Y,Matrix):
        M:=Matrix(Y):
        # lprint("Matrix M: ",M):
        row,col:=Dimension(M):
        lprint("row: ",row, " col: ",col):
        
        if row =1 then 
            lin_sys:=false:
            u:=Interp(alpha,Y,x)mod p:
            # lprint("NDSA: Single equation case - u: ",u):
            result:=[[MQRFR(m,u,0,1,p)]]:
            # lprint("result =",result):
            dq:=result[1][3]:
            print("NDSA: dq: ",dq):
            # print("f=",result[1][1]):
            # print("g=",result[1][2]):

        else
            lin_sys:=true: 
            u:=get_u(M,col,alpha,p):
        end if:
        # lprint("NDSA: u: ",u):
        # print("nops u: ",nops(u));  
        
        if lin_sys = true then  
            for i from 1 to nops(u) do 
            # put a check for MQRFR_done
                if(MQRFR_done[i]) then 
                    print("NDSA: Skipping MQRFR for equation ",i," as already done"):
                    next:
                end if;
                temp:=[op(temp),MQRFR(m,u[i],0,1,p)]:
                # print("NDSA: result temp before assignment: ",result):
                # print("NDSA: temp",i,": ",temp):
                # result:=[op(result),temp]:
                result[i]:=temp:
                print("NDSA: result",i,": ",result[i]):
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
            print("NDSA: mqrfr_status: ",MQRFR_done):
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
    if(T>2^5)then break; end if; # Safety break to avoid infinite loops during testing
    end do:
end proc:


