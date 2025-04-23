
generate_powers:=proc(t_old,T,point_,num_var,p)
    local i,j:
    print("In generate_powers");
    print("T=",T);
    print("point=",point_);
    return modp([seq([seq(point_[j]^i,j = 1..num_var)], i = t_old..2*T)],p):
end proc:  

get_lambda:=proc(B,anchor_point,shift_,T,num_var,p)# Needs correction find a way to keep roots in mem for t 
# and t+1 to check if the number of roots is equal to the number of terms. 
    print("In get_lambda");
    local t,flag,Y,Lambda,terms,R,num_points:
    # f:=f*c mod p:
    # g:=g*c mod p:
    # Note: t_old =1 to 2*t was causing an error in Zippel VanderMonde Solver. 
    t_old:=0;
    t:=T:
    print("t=",t):
    # anchor_points:=generate_powers(t_old,t,anchor_point,num_var,p):
    memo:=[]:
    i:=1;
    num_points:=t:
    berlekamp_failure:=true:
    numerator_failure:=true:
    denominator_failure:=true:
    # counter:=0;
    num_vals:=[];
    den_vals:=[];
    j_init:=2:
    while true do

        # terms_num:=[];
        # terms_den:=[]; 
        print("i=",i):
        print("t=",t):
        print("t_old=",t_old):  
        print("anchor_points: ",anchor_point):
        # [2^j,3^j]
        # for ii from t_old to 2*t do 
        for ii from t_old to 2*t-1 do 
           memo:=[op(memo),[seq(anchor_point[j]^ii,j=1..num_var)]] mod p:
        end do;
        print("memo=",memo):
        # break;
        # anchor_points:=[op(anchor_points),[seq(modp(anchor_point[ii]^i,p),ii=1..num_var)]]:
        print("numelems(memo) = ",numelems(memo));
        # get degree here 
        if berlekamp_failure and i=1 then
            print("IN FIRST ITERATION"): 
            f,g,lg,num_points:=evaluation_num_den(B,memo[1],shift_,num_var,p,num_points): 
            c:=1/eval(f,x=0) mod p;
            num_eval:=eval(f,x=memo[1][1])*c mod p:
            num_vals:=[op(num_vals),num_eval]:
            den_eval:=eval(g,x=memo[1][1])*c mod p:
            den_vals:=[op(den_vals),den_eval]:
            deg_num:=degree(f,x):
            deg_den:=degree(g,x):
            print("f= ",f):
            print("g= ",g):
            num_points:=deg_num+deg_den+2:
            # For this to work dq in 3_eval.mpl needs to be > 1 on line 40
            print("num_points= ",num_points): 
        end if:
        print("j_init=",j_init):
    
        for j from j_init to numelems(memo) do 
            print("Back in get_lambda j=",j):
            # f,g,lg:=evaluation_num_den(B,memo[j][1],shift_,num_var,p):
            print("memo[j]=",memo[j]):
            f,g,lg:=evaluation_num_den(B,memo[j],shift_,num_var,p,num_points): 
            c:=1/eval(f,x=0) mod p;
            print("f= ",f):
            print("g= ",g):
            if numerator_failure then 
                # print("numerator_failure= ",numerator_failure):
                num_eval:=eval(f,x=memo[j][1])*c mod p:
                print("num_eval=",num_eval):
                num_vals:=[op(num_vals),num_eval]:
                print("num_vals= ",num_vals):
            end if:
            # true_num:=[op(true_num),eval(ff,{x=memo[j][1],y=memo[j][2]}) mod p]:
            # print("true_num= ",true_num):
            if denominator_failure then 
                # print("denominator_failure= ",denominator_failure):
                den_eval:=eval(g,x=memo[j][1])*c mod p:
                # print("den_eval=",den_eval):
                den_vals:=[op(den_vals),den_eval]:               
                print("den_vals= ",den_vals):
            end if:
            print("________________________________________________________________"):
        end do:

        # break;
        print("i=",i):
        # Call BMEA only on the polynomial which caused the failure, while keeping the value of the other polynomial.
        if numerator_failure then 
            lambda_num := BMEA(num_vals,p,Z):
            terms_num:=degree(lambda_num,Z):
            # terms_num[i]:=degree(lambda_num,Z):
            R_num := Roots(lambda_num) mod p:
            print("lambda_num= ",lambda_num):
             print("terms_num[i]= ",terms_num[i]):
            print("R_num= ",R_num):
        end if:
        if denominator_failure then 
            lambda_den := BMEA(den_vals,p,Z):
            print("lambda_den= ",lambda_den): 
            terms_den:=degree(lambda_den,Z):
            # terms_den[i]:=degree(lambda_den,Z):
            print("terms_den[i]= ",terms_den[i]):
            # R_num := Roots(lambda_num) mod p:
            R_den := Roots(lambda_den) mod p:
        end if:
       
          
        # terms_num[i]:=degree(lambda_num,Z):
        # terms_den[i]:=degree(lambda_den,Z):
        # print("terms_den[i]= ",terms_den[i]):
        # R_num := Roots(lambda_num) mod p:
        # R_den := Roots(lambda_den) mod p:
        # print("R_num= ",R_num):
        # print("R_den= ",R_den):

        print("num of roots of lambda_num =",nops(R_num)):
        print("num of roots of lambda_den =",nops(R_den)):
        # print("i=",i):
        # break:
        #  Find out if failure is caused due to the numerator or the denominator.
        # if nops(R_num)=terms_num[i] then 
        if nops(R_num)=terms_num then 
            print("Numerator Success!!"):
            numerator_failure:=false:
        end if:
        # if nops(R_den)=terms_den[i] then 
        if nops(R_den)=terms_den then 
            print("Denominator Success!!"):
            denominator_failure:=false:
        end if:
        # if i=1 and terms_num[i] = nops(R_num) and terms_den[i] = nops(R_den) then # Success condition 1
        #     return terms_num[i],terms_den[i],lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:
        # # elif i > 1 and (nops(R_num)<terms_num[i] or nops(R_den)<terms_den[i]) then 

        # # else Break condition needs work. 
        # elif i>1 and terms_num[i-1]=terms_num[i] and terms_num[i] = nops(R_num) # success condition 2
        # and terms_den[i-1]=terms_den[i] and terms_den[i] = nops(R_den) then 
        #     # print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");
        #         return terms_num[i],terms_den[i],lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:
        # end if: 
        # end if:
        print("terms_num[i]= ",terms_num):
        print("terms_den[i]= ",terms_den):
        berlekamp_failure:=numerator_failure or denominator_failure:
        print("berlekamp_failure= ",berlekamp_failure):
        if not(berlekamp_failure) then 
            print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");  
            return terms_num,terms_den,lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:
        end if:
        i:=i+1:
        t_old:=2*t:
        j_init:=t_old+1:
        t:=t*2:
        
        print("t=",t):
        print("t_old=",t_old):

        print("Berlekamp Massey failure!!"):
        
        if(i = 5 )then 
            break;
        end if:
        break;
    end do:
    # print("terms=",terms[i]):
end proc: