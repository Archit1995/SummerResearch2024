with(NumberTheory): 
generate_powers:=proc(t_old,T,point_,num_var,p)
    local i,j:
    print("In generate_powers");
    print("T=",T);
    print("point=",point_);
    return modp([seq([seq(point_[j]^i,j = 1..num_var)], i = t_old..2*T)],p):
end proc:  

get_lambda:=proc(B,anchor_point,shift_,T,num_var,p,test_case)# Needs correction find a way to keep roots in mem for t 
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
    numerator_done:=false:
    denominator_done:=false:
    # counter:=0;
    num_vals:=[];
    den_vals:=[];
    true_num:=[];
    true_den:=[];
    num_ratios:=[];
    den_ratios:=[];
    raw_num:=[];
    raw_den:=[];
    j_init:=2:

    f,g,lg,num_points:=evaluation_num_den(B,[seq(1,i=1..num_var)],shift_,num_var,p,num_points): 
    num_eval:=eval(f,x=1):
    num_vals:=[op(num_vals),num_eval]:
    den_eval:=eval(g,x=1):
    den_vals:=[op(den_vals),den_eval]:
    deg_num:=degree(f,x):
    deg_den:=degree(g,x):
    num_points:=deg_num+deg_den+2:
    print("num_points= ",num_points):

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
        
        print("j_init=",j_init):
    
        for j from j_init to numelems(memo) do 
            print("Back in get_lambda j=",j):
            # f,g,lg:=evaluation_num_den(B,memo[j][1],shift_,num_var,p):
            print("memo[j]=",memo[j]):
            f,g,lg:=evaluation_num_den(B,memo[j],shift_,num_var,p,num_points): 
            # Ratrecon 
            # m,u,f,g,lg:=evaluation_num_den(B,memo[1],shift_,num_var,p,num_points): 
            c:=1/eval(f,x=0) mod p;
            # ratrecon_result:=Ratrecon(u, m, x,deg_num,deg_den, 'nn', 'dd') mod p:
            # print("ratrecon result = ",ratrecon_result):
            # print("nn = ", nn);
            # print("dd = ",dd);
            # print("-------------------------------------------------------------");
            print("-----------------------------------------------------"):
            # y= shift_[1]*(x-memo[1][1])+memo[1][2]:
            true_f:=simplify(subs(y=shift_[1]*(x-memo[j][1])+memo[j][2],ff))mod p:
            true_g:=simplify(subs(y=shift_[1]*(x-memo[j][1])+memo[j][2],gg))mod p:
            # print("true f = ",true_f):
            # print("true g = ",true_g):
            # print("true num value= ",eval(true_f,x=memo[j][1]) mod p):
            # print("true den value= ",eval(true_g,x=memo[j][1]) mod p):

            # print("--------------------------------------"):
            # print("c = 1/f(0) = ",c):
            # print("f*c mod p= ",f*c mod p):
            # print("g*c mod p= ",g*c mod p):
            # print("f= ",f):
            # print("g= ",g):
            # print("gcd = ",Gcd(f,g) mod p);
            # _resultant:=Resultant(f,g,x) mod p:
            # print("Resultant= ",_resultant):
            # print("lg= ",lg):
            # print("MultiplicativeOrder(lg,p)=",MultiplicativeOrder(lg,p)):
            # print("c= ",c):
            # print("lg/ratio_= ",lg/num_ratios[1] mod p):
            # print("ratio/lg= ",num_ratios[1]/lg mod p):
            # print("shift/lg= ",shift_[1]/lg mod p):
            # print("shift*lg= ",shift_[1]*lg mod p):


            if not(numerator_done) then 
                # print("numerator_done= ",numerator_done):
                num_eval:=eval(f,x=memo[j][1]):
                # print("num_eval=",num_eval):
                num_vals:=[op(num_vals),num_eval]:
                # print("num_vals= ",num_vals):
                raw_num:=[op(raw_num),eval(f,x=memo[j][1]) mod p]:
                true_num:=[op(true_num),eval(ff,{x=memo[j][1],y=memo[j][2]}) mod p]:
                # print("raw_num= ",raw_num):
                # print("true_num= ",true_num):
                # num_ratios:=[op(num_ratios),raw_num[j]/true_num[j] mod p]:
                # print("num_ratios= ",num_ratios):
                # print("f(0) = ",1/c mod p);
                # print("f(1) = ",eval(f,x=1) mod p);
            end if:

            
            if not(denominator_done) then 
                # print("denominator_done= ",denominator_done):
                den_eval:=eval(g,x=memo[j][1]):
                # print("den_eval=",den_eval):
                den_vals:=[op(den_vals),den_eval]:               
                print("den_vals= ",den_vals):
                # true_den:=[op(true_den),eval(gg,{x=memo[j][1],y=memo[j][2]}) mod p]:
                # raw_den:=[op(raw_den),eval(g,x=memo[j][1]) mod p]:
                # print("raw_den= ",raw_den):
                # print("true_den= ",true_den):
                # den_ratios:=[op(den_ratios),raw_den[j]/true_den[j] mod p]: 
                # print("den_ratios= ",den_ratios):
                # print("f(0) = ",1/c mod p);
            end if:

            print("________________________________________________________________"):
        end do:

        # break;
        print("i=",i):
        # Call BMEA only on the polynomial which caused the done, while keeping the value of the other polynomial.
        if not(numerator_done) then 
            lambda_num := BMEA(num_vals,p,Z):
            terms_num:=degree(lambda_num,Z):
            # terms_num[i]:=degree(lambda_num,Z):
            R_num := Roots(lambda_num) mod p:
            print("lambda_num= ",lambda_num):
            print("terms_num[i]= ",terms_num):
            print("R_num= ",R_num):
        end if:
        if not(denominator_done) then 
            lambda_den := BMEA(den_vals,p,Z):
            print("lambda_den= ",lambda_den): 
            terms_den:=degree(lambda_den,Z):
            # terms_den[i]:=degree(lambda_den,Z):
            print("terms_den[i]= ",terms_den):
            # R_num := Roots(lambda_num) mod p:
            R_den := Roots(lambda_den) mod p:
            print("R_den= ",R_den):
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
        #  Find out if done is caused due to the numerator or the denominator.
        # if nops(R_num)=terms_num[i] then 
        if nops(R_num)=terms_num then 
            print("Numerator Success!!"):
            numerator_done:=true:
        end if:
        # if nops(R_den)=terms_den[i] then 
        if nops(R_den)=terms_den then 
            print("Denominator Success!!"):
            denominator_done:=true:
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
        if numerator_done and denominator_done then 
            print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");  
            # print("terms_num= ",terms_num):
            # print("terms_den= ",terms_den):
            # print("lambda_num= ",lambda_num):
            # print("lambda_den= ",lambda_den):
            # print("R_num= ",R_num):
            # print("R_den= ",R_den):
            # print("num_vals= ",num_vals):
            # print("den_vals= ",den_vals):
            return terms_num,terms_den,lambda_num,lambda_den,R_num,R_den,num_vals,den_vals:
        end if:
        i:=i+1:
        t_old:=2*t:
        j_init:=t_old+1:
        t:=t*2:
        
        print("t=",t):
        print("t_old=",t_old):

        print("Berlekamp Massey done!!"):
        
        # if test_case = 3 or test_case = 11 then 
        #     if(i = 3 )then 
        #         return terms_num,terms_den,lambda_num,lambda_den,R_num,R_den,num_vals,den_vals,den_ratios[1]:
        #         # break;
        #     end if:
        # end if:
    end do:
    # print("terms=",terms[i]):
end proc: