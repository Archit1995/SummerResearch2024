read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
read "./4_BMEA.mpl":
p:=2^31-1:
T:=4:
# p:=811:
# p:=19:
ff:=x*y+1:
gg:=x+y:
# ff:=x^3*y^3+x^2*y^2+x*y+1:
# gg:=x^2+y^2+x+y+3:
num_var:=2:
vars:={x,y}:
test:=8:
B:=Construct_Rational_Blackbox(ff,gg,vars):
anchor_point:=[2,3]:
shift_:=[10]:
# j:=1:
for j from 1 to 10 do 
    print("j= ",j):
    anchor_points:=[seq(modp(anchor_point[i]^j,p),i=1..num_var)]:
    print("anchor_points= ",anchor_points):
    f,g,lg:=evaluation_num_den(B,anchor_points,shift_,num_var,p):
    print("f= ",f):
    print("g= ",g):
    c:=1/eval(f,x=0) mod p;
    print("anchor_points[1] = ",anchor_points[1]):
    num_eval[j]:=eval(f,x=anchor_points[1])*c mod p:
    den_eval[j]:=eval(g,x=anchor_points[1])*c mod p:
    num_vals:=convert(num_eval,list):
    print("num_vals= ",num_vals):
    true_num[j]:=eval(ff,{x=anchor_points[1],y=anchor_points[2]}) mod p:
    print("true_num= ",true_num):
    den_vals:=convert(den_eval,list):
    true_den[j]:=eval(gg,{x=anchor_points[1],y=anchor_points[2]}) mod p:
    print("true_den= ",true_den):
    print("den_vals= ",den_vals):
    print("=====================================================");
    lambda_num := BMEA(num_vals,p,Z):
    lambda_den := BMEA(den_vals,p,Z):
    terms_num[j]:=degree(lambda_num,Z):
    terms_den[j]:=degree(lambda_den,Z):
    R_num := Roots(lambda_num) mod p:
    R_den := Roots(lambda_den) mod p:
    break;
    if nops(R_num)<terms_num[j] and nops(R_den)<terms_den[j] then 
        j:=j*2:
    elif terms_num[j-1]=terms_num[j] and terms_num[j] = nops(R_num) and terms_den[j-1]=terms_den[j] and terms_den[j] = nops(R_den) then 
        break:
    end if:
    # print("====================================================="):
    # print("Estimated c/ acutal c = ",lg/actual_c mod p);
    # print("Estimated 1/c/ acutal 1/c = ",1/lg/actual_c_inv mod p);
    # print("Estimated 1/c - acutal 1/c = ",1/lg-actual_c_inv mod p);
    # print("___________________________________________________________________");
end do:
# print("lambda_num= ",lambda_num):
# print("lambda_den= ",lambda_den):
