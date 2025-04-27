read "./0data_gen.mpl":
read "./1_construct_BB.mpl":
read "./2_mqrfr.mpl":
read "./3_eval.mpl":
with(NumberTheory):
vars:={x}:
# p:=811:
p:=2^31-1:
deg:=[7,3];
# p:=31:
T:=4:




# print("MultiplicativeOrder of actual_lc", MultiplicativeOrder(actual_lc,p)):
# print("MultiplicativeOrder of mqrfr_lc", MultiplicativeOrder(mqrfr_lc,p)):
mqrfr_lcoeff:=[]:
# ff,gg,actual_lc:=data_generator(1,deg,p):
for i from 1 to 10 do
    print("i= ",i):
    ff,gg,actual_lc:=data_generator(1,deg,p):
    print("ff= ",ff):
    print("gg= ",gg):
    B:=Construct_Rational_Blackbox(ff,gg,vars):
    f,g,mqrfr_lc:=evaluation_num_den(B,T,p):
    mqrfr_lcoeff:=[op(mqrfr_lcoeff),mqrfr_lc]:
    print("f= ",f):
    print("g= ",g):
    print("actual_lc= ",actual_lc):
    print("mqrfr_lc= ",mqrfr_lc):
    print("MultiplicativeOrder of actual_lc", MultiplicativeOrder(actual_lc,p)):
    print("MultiplicativeOrder of mqrfr_lc", MultiplicativeOrder(mqrfr_lc,p)):
    
    print("___________________________________________________________________"):
end do:
mqrfr_lcoeff;
# print("mqrfr_lcoeff= ",mqrfr_lcoeff):
# with(plots):
# listplot(mqrfr_lcoeff):
# actual_lc_inv:=1/actual_lc mod p:
# mqrfr_lcoeff*actual_lc_inv mod p;

# evaluation_num_den(B,T,p):

# print("lcg/lcggg= ",lcg/lcggg mod p):
# print("lcggg/lcg= ",lcggg/lcg mod p):