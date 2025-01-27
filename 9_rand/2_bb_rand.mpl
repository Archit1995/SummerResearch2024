Construct_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        local u,v,var,num,den,a;
        var:=vars:
        num:=f:
        denom_:=g:
        a:=num/denom_:
        # print("a: ",a);
        # return [seq(Eval(a,x=i) mod p,i=point_[1..numelems(point_)])]:
        return Eval(a,x=point_) mod p:
    end proc:
    return BB:
end proc:
# generate_random_vector:=proc(n,p)
#     r:=rand(p);
#     return [seq(r(),i=1..n)];
# end proc:
p:=2^31-1;
r:=rand(p);
num_deg:=17:
denom_deg:=15:
f0:=randpoly(x,degree=num_deg):
#  lc_f0:=lcoeff(f0):
#  f0:=f0/lc_f0 mod p:
f1:=randpoly(x,degree=denom_deg):
lc_f1:=lcoeff(f1):
f1:=f1/lc_f1 mod p:
vars:={x}:
B:=Construct_Blackbox(f0,f1,vars):
[seq(B(r(),p),i=1..10)];