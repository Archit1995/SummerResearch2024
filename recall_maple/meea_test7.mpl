B:=proc(a,b,point_)
    local c;
    c:= a/b;
    return [seq(eval(c,x=i),i=point_[1..numelems(point_)])];
end proc;
eea_gcd:=proc(r0,r1,s0,t0,s1,t1)
   local r2,q,s2,t2;
    if r1 = 0 then return r0, s0, t0 end if;
    # 1. find quotient and remainder
    q:= quo(r0,r1,x,'r2');
    print("q=",q);
    print("r2=",r2);
    #  reduction step to update s and t.
    s2:=s0-q*s1;
    t2:=t0-q*t1:
    # Normalization
    if r2<> 0 then 
        rho:=lcoeff(r2);
        r2:=r2/rho;
        s2:=s2/rho;
        t2:=t2/rho;
    end if;
    eea_gcd(r1,r2,s1,t1,s2,t2);
end proc;
num_deg:=20;
denom_deg:=10;
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

g,s,t:=eea_gcd(f0/rho0,f1/rho1,1/rho0,0,0,1/rho1);
g-expand(s*f0+t*f1);