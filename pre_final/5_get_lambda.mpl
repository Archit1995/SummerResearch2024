
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
    local t,flag,Y,Lambda,terms,R:
    # f:=f*c mod p:
    # g:=g*c mod p:
    # Note: t_old =1 to 2*t was causing an error in Zippel VanderMonde Solver. 
    t_old:=0;
    t:=T:
    # anchor_points:=generate_powers(t_old,t,anchor_point,num_var,p):
    memo:=[]:
    i:=1;
    # counter:=0;
    while true do
        num_vals:=[];
        den_vals:=[];
        # terms_num:=[];
        # terms_den:=[]; 
        print("i=",i):
        print("t=",t):
        print("t_old=",t_old):  
        print("anchor_points: ",anchor_points):
        # [2^j,3^j]
        # for ii from t_old to 2*t do 
        for ii from t_old to 2*t-1 do 
           memo:=[op(memo),[seq(anchor_point[j]^ii,j=1..2)]] mod p:
        end do;
        print("memo=",memo):
        # break;
        # anchor_points:=[op(anchor_points),[seq(modp(anchor_point[ii]^i,p),ii=1..num_var)]]:
        print("numelems(memo) = ",numelems(memo));
        for j from 1 to numelems(memo) do 
            print("Back in get_lambda j=",j):
            # f,g,lg:=evaluation_num_den(B,memo[j][1],shift_,num_var,p):
            f,g,lg:=evaluation_num_den(B,memo[j],shift_,num_var,p):  
            c:=1/eval(f,x=0) mod p;
            print("f= ",f):
            print("g= ",g):
            print("=memo[j][1]=",memo[j][1]):
            num_eval:=eval(f,x=memo[j][1])*c mod p:
            print("num_eval=",num_eval):
            num_vals:=[op(num_vals),num_eval]:
            print("num_vals= ",num_vals):
            # true_num:=[op(true_num),eval(ff,{x=memo[j][1],y=memo[j][2]}) mod p]:
            # print("true_num= ",true_num):
            den_eval:=eval(g,x=memo[j][1])*c mod p:
            # print("den_eval=",den_eval):
            den_vals:=[op(den_vals),den_eval]:
            # true_den:=[op(true_den),eval(gg,{x=memo[j][1],y=memo[j][2]}) mod p]:
            # print("true_den= ",true_den):
            print("den_vals= ",den_vals):
            print("________________________________________________________________"):
        end do:

        # break;
        print("i=",i):
        lambda_num := BMEA(num_vals,p,Z):
        lambda_den := BMEA(den_vals,p,Z):
        print("lambda_num= ",lambda_num):
        print("lambda_den= ",lambda_den):   
        terms_num[i]:=degree(lambda_num,Z):
        terms_den[i]:=degree(lambda_den,Z):
        print("terms_num[i]= ",terms_num[i]):
        print("terms_den[i]= ",terms_den[i]):
        R_num := Roots(lambda_num) mod p:
        R_den := Roots(lambda_den) mod p:
        print("R_num= ",R_num):
        print("R_den= ",R_den):

        print("num of roots of lambda_num =",nops(R_num)):
        print("num of roots of lambda_den =",nops(R_den)):
        # print("i=",i):
        # break:
        if i=1 and terms_num[i] = nops(R_num) and terms_den[i] = nops(R_den) then # Success condition 1
            return terms_num[i],terms_den[i],lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:
        # elif i > 1 and (nops(R_num)<terms_num[i] or nops(R_den)<terms_den[i]) then 

        # else Break condition needs work. 
        elif i>1 and terms_num[i-1]=terms_num[i] and terms_num[i] = nops(R_num) # success condition 2
        and terms_den[i-1]=terms_den[i] and terms_den[i] = nops(R_den) then 
            # print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");
                return terms_num[i],terms_den[i],lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:
        end if: 
        # end if:
        i:=i+1:
        t_old:=2*t+1:
        t:=t*2:
        print("Berlekamp Massey failure!!"):
        if(i = 5 )then 
            break;
        end if:
        # break;
    end do:
    # print("terms=",terms[i]):
end proc: