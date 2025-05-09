read "./0_data_gen.mpl":
read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
read "./4_BMEA.mpl":
num_var,vars,p,T,ff,gg:=data_generator(3):
test:=8:
B:=Construct_Rational_Blackbox(ff,gg,vars):
anchor_point:=[2,3]:
shift_:=[10]:
t_old:=0;
t:=T:
memo:=[]:
i:=1;
while true do
    num_vals:=[];
    den_vals:=[];
    true_num:=[];
    true_den:=[]; 
    print("i=",i):
    print("t=",t):
    print("t_old=",t_old):  
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
        f,g,lg:=evaluation_num_den(B,memo[j],shift_,num_var,p):  
        c:=1/eval(f,x=0) mod p;
        print("f= ",f):
        print("g= ",g):
        num_eval:=eval(f,x=memo[j][1])*c mod p:
        # print("num_eval=",num_eval):
        num_vals:=[op(num_vals),num_eval]:
        print("num_vals= ",num_vals):
        true_num:=[op(true_num),eval(ff,{x=memo[j][1],y=memo[j][2]}) mod p]:
        print("true_num= ",true_num):
        den_eval:=eval(g,x=memo[j][1])*c mod p:
        # print("den_eval=",den_eval):
        den_vals:=[op(den_vals),den_eval]:
        true_den:=[op(true_den),eval(gg,{x=memo[j][1],y=memo[j][2]}) mod p]:
        print("true_den= ",true_den):
        print("den_vals= ",den_vals):
        print("________________________________________________________________"):
    end do:    
    i:=i+1:
    t_old:=2*t+1:
    t:=t*2:
    if(i = 5 )then 
        break;
    end if:
    # break;
end do:


# for j from 1 to 10 do 
#     print("j= ",j):
#     anchor_points:=[seq(modp(anchor_point[i]^j,p),i=1..num_var)]:
#     print("anchor_points= ",anchor_points):
#     f,g,lg:=evaluation_num_den(B,anchor_points,shift_,num_var,p):
#     print("f= ",f):
#     print("g= ",g):
#     c:=1/eval(f,x=0) mod p;
#     print("anchor_points[1] = ",anchor_points[1]):
#     num_eval[j]:=eval(f,x=anchor_points[1])*c mod p:
#     den_eval[j]:=eval(g,x=anchor_points[1])*c mod p:
#     num_vals:=convert(num_eval,list):
#     print("num_vals= ",num_vals):
#     true_num[j]:=eval(ff,{x=anchor_points[1],y=anchor_points[2]}) mod p:
#     print("true_num= ",true_num):
#     den_vals:=convert(den_eval,list):
#     true_den[j]:=eval(gg,{x=anchor_points[1],y=anchor_points[2]}) mod p:
#     print("true_den= ",true_den):
#     print("den_vals= ",den_vals):
#     print("=====================================================");
#     lambda_num := BMEA(num_vals,p,Z):
#     lambda_den := BMEA(den_vals,p,Z):
#     terms_num[j]:=degree(lambda_num,Z):
#     terms_den[j]:=degree(lambda_den,Z):
#     R_num := Roots(lambda_num) mod p:
#     R_den := Roots(lambda_den) mod p:
#     break;
#     if nops(R_num)<terms_num[j] and nops(R_den)<terms_den[j] then 
#         j:=j*2:
#     elif terms_num[j-1]=terms_num[j] and terms_num[j] = nops(R_num) and terms_den[j-1]=terms_den[j] and terms_den[j] = nops(R_den) then 
#         break:
#     end if:
# end do:

