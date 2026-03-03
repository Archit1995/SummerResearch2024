########################################################
# 2. NDSA - For Vector Black Boxes
########################################################

with(LinearAlgebra):

Deterministic_NDSA:=proc(B,sigma_,beta_,num_var,p,num_points,max_points,num_eqn)
    local correct_degree,T,alpha,m,phi_,Psi_alpha,Y,u,f,g,dq,lcg,np,nv,i,r,
          lin_sys,temp,result,count,M,row,col,DQ,MQRFR_done:
    print("In Deterministic_NDSA");
    # print("Deterministic_NDSA: num_eqn:= ",num_eqn);
    # print("Deterministic_NDSA: num_points:= ",num_points);
    lin_sys:=false:
    T:=max_points:
    temp:=[]:
    # result:=[]:
    result:=[seq([], i=1..num_eqn)]:
    # result:=Array(1..num_eqn,[]):
    count:=0:
    r:=rand(p):
    # alpha:=[seq(r()+r() mod p,i=1..T)]:
    alpha:=[seq(i ,i=2..T+1)]:
    print("Deterministic_NDSA: alpha: ",alpha):
    Psi_alpha:=get_point_on_affine_line(num_var,alpha,beta_,sigma_,p,T):
    print("Deterministic_NDSA:Psi_alpha: ",Psi_alpha):
    Y:=[seq(B(Psi_alpha[i],p),i=1..T)]:
    print("Deterministic_NDSA: Y: ",Y);
    for i from 1 to num_eqn do
        # print("Deterministic_NDSA: i= ",i,"num_points[",i,"]= ",num_points[i]); 
        # print(",alpha[1..num_points[i]]= ",alpha[1..num_points[i]]);
        m:=expand(product(x-alpha[j],j=1..num_points[i])) mod p:
        print("Deterministic_NDSA:m: ",m):
        M:=convert(Y[1..num_points[i]],Matrix):
        row,col:=Dimension(M):
        # print("row: ",row, " col: ",col):
        # print("Deterministic_NDSA: Matrix M: ",M):
        if row =1 then 
            u:=Interp(alpha[1..num_points[i]],Y[1..num_points[i]],x)mod p:
        else
            lin_sys:=true: 
            # u:=get_u(M,col,alpha[1..num_points[i]],p):
            u:=Deterministic_get_u(M,i,alpha[1..num_points[i]],p):
        end if:
        lprint("Deterministic_NDSA: u: ",u):    
            if lin_sys = false then  
            result:=[MQRFR(m,u,0,1,p)]:
            dq:=result[1][3]:
        else 
            result[i]:=[MQRFR(m,u,0,1,p)]:
            # lprint("Deterministic_NDSA: result",i,": ",result[i]):
            # for ii from 1 to nops(u) do 
            #     result[ii]:=[MQRFR(m,u[ii],0,1,p)]:
            #     print("Deterministic_NDSA: result",ii,": ",result[ii]):
            # end do:
        end if;
        # print("_____________________________________________________________"):
    end do:
    # lprint("Deterministic_NDSA:entire result: ",result):
    return result:
end proc:


