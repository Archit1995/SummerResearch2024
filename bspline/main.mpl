read "./0_constructBB.mpl":
read "./reordering.mpl":
read "./data_gen.mpl":
read "./Kaltofen_Yang_header.mpl":
Sys,Vars,params,p:=get_data("bspline"):
nvar:=nops(Vars);
nparam:=nops(params);
LBB:=Constuct_Sys_Blackbox(Sys,Vars,params):
# MRFI(LBB,nparam,params,p):
sigma_:=[]:
Primes:=[seq(ithprime(i),i=1..nparam)]:
for j from 1 to 4 do 
    sigma_:=[op(sigma_),[seq(Primes[i]^j mod p,i=1..nops(Primes))]]:
    print("sigma_: ",sigma_[j]);
    print("LBB: ",LBB(sigma_[j],p)):
    print(op(4,eval(LBB)));
    print("__________________________________"):
end do:
# sigma_;   
    # LBB(sigma_,p);
# eval(LBB);
# op(4,eval(LBB))[sigma_[j],p][i] =xi=fi(y1..ym)/gi(y1..ym)
# op(4,eval(LBB))[sigma_[1],p][1];
# for i from 1 to nops(eval(LBB)) do 
#     print("i =",i,"op(i) = ",);
# end do;
