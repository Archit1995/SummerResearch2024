read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
p:=2^31-1:
# p:=811:
# p:=19:
ff:=x*y+1:
gg:=x+y:
num_var:=2:
vars:={x,y}:
test:=8:
B:=Construct_Rational_Blackbox(ff,gg,vars):
anchor_point:=[2,3]:
shift_:=[5]:
for j from 1 to test do
    print("j= ",j):
    anchor_points:=[seq(modp(anchor_point[i]^j,p),i=1..2)]:
    print("anchor_points= ",anchor_points):
    f,g,lg:=evaluation_num_den(B,anchor_points,shift_,num_var,p):
    print("f= ",f):
    print("g= ",g):

    f_2:=eval(f,x=anchor_points[1])mod p:
    g_2:=eval(g,x=anchor_points[1])mod p:
    # print("f_2^",j,"= ",f_2):
    # print("g_2^",j,"= ",g_2):
    ff_2_3:=eval(ff,{x=anchor_points[1],y=anchor_points[2]})  mod p:
    gg_2_3:=eval(gg,{x=anchor_points[1],y=anchor_points[2]})  mod p:
    # print("ff_2_3^",j,"= ",ff_2_3):
    # print("gg_2_3^",j,"= ",gg_2_3):
    actual_c:=ff_2_3/f_2 mod p:
    actual_c_inv:=f_2/ff_2_3 mod p:
    print("====================================================="):
    print("calculated/estimated c= ",lg):
    # print("actual c = ",ff_2_3/f_2 mod p);
    print("actual c = ",gg_2_3/g_2 mod p);
    print("calculated/estimated 1/c= ",1/lg mod p);
    # print("actual 1/c = ",f_2/ff_2_3 mod p);
    print("actual 1/c = ",g_2/gg_2_3 mod p);

    print("====================================================="):
    print("Estimated c/ acutal c = ",lg/actual_c mod p);
    print("Estimated 1/c/ acutal 1/c = ",1/lg/actual_c_inv mod p);
    print("Estimated 1/c - acutal 1/c = ",1/lg-actual_c_inv mod p);
    print("___________________________________________________________________");
end do: