B:=proc(a,b,point_)
    local c;
    c:= a/b;
    return [seq(eval(c,x=i),i=point_[1..numelems(point_)])];
end proc;
eea_gcd:=proc(r0,r1,s0,t0,s1,t1,p)
   local r2,q,s2,t2;
    if r1 = 0 then return r0, s0, t0 end if;
    # 1. find quotient and remainder
    q:= Quo(r0,r1,x,'r2') mod p;
    print("q=",q);
    print("r2=",r2);
    #  reduction step to update s and t.
    s2:=s0-q*s1 mod p;
    t2:=t0-q*t1 mod p;:
    # Normalization
    if r2<> 0 then 
        rho:=lcoeff(r2);
        r2:=r2/rho mod p;
        s2:=s2/rho  mod p;  
        t2:=t2/rho mod p;
    end if;
    eea_gcd(r1,r2,s1,t1,s2,t2,p);
end proc;
num_deg:=2;
denom_deg:=1;
# d:=num_deg+denom_deg;
# num_points:=d+1;
# deg m = d > N+M
# 1. generate the points
# r:=rand(1000..10000);
# points:=[seq(r(),1..num_points)];
# m:expand(product(x-points[i],i=1..num_points));
# 2. generate the polynomial
f0:=randpoly(x,degree=num_deg);
rho0:=lcoeff(f0);
f1:=randpoly(x,degree=denom_deg);
rho1:=lcoeff(f1);
# y:=B(f1,f2,points);
# u:=interp(y,points,x);
p=:19;
rho0_inv:=1/rho0 mod p;
rho1_inv:=1/rho1 mod p;
g,s,t:=eea_gcd(f0*rho0_inv,f1*rho1_inv,rho0_inv,0,0,rho1_inv,p);
g-Expand(s*f0+t*f1) mod p;