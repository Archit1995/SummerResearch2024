B:=proc(a,b,p,point_)
    local c;
    c:= a/b;
    return [seq(Eval(c,x=i) mod p,i=point_[1..numelems(point_)])];
end proc;
eea_gcd:=proc(r0,r1,s0,t0,s1,t1,p)
   local r2,q,s2,t2;
    if r1 = 0 then return r0, s0, t0 end if;
    # 1. find quotient and remainder
    q:= Quo(r0,r1,x,'r2') mod p;
    print("q=",q);
    # print("r2=",r2);
    #  reduction step to update s and t.
    # s2:=s0-q*s1 mod p;
    t2:=t0-q*t1 mod p;
       
    # Normalization
    if r2<> 0 then 
        rho:=lcoeff(r2);
        r2:=r2/rho mod p;
        print("r2=",Expand(r2) mod p);
        # s2:=s2/rho mod p;  
        t2:=t2/rho mod p;
        print("t2=",Expand(t2)mod p); 
        print("--------------------");
    end if;

    eea_gcd(r1,r2,s1,t1,s2,t2,p);
end proc;
p:=19;
# num_deg:=5;
# denom_deg:=3;
# d:=num_deg+denom_deg;
# num_points:=d+1;
# # deg m = d > N+M
# # 1. generate the points
# r:=rand(1000..10000);
# points:=[seq(r(),1..num_points)];
# m:=expand(product(x-points[i],i=1..num_points));
# lc_m:=lcoeff(m);
# m_m:=m/lc_m mod p;

# 2. generate the polynomial
# f0:=randpoly(x,degree=num_deg);
# rho3:=lcoeff(f0);
# f0:=f0/rho3 mod p;
# f1:=randpoly(x,degree=denom_deg);
# rho4:=lcoeff(f1);
# f1:=f1/rho4 mod p;
f0:=3*x^2+4;
f1:=x+2;
points:=[15,13,6,2];
y:=B(f0,f1,p,points);
u:=Interp(y,points,x) mod p;
lc_u:=lcoeff(u);
m_u:=u/lc_u mod p;

m:=expand(product(x-points[i],i=1..numelems(points)));
lc_m:=lcoeff(m);
m_m:=m/lc_m mod p;

lc_m_inv:=1/lc_m mod p;
lc_u_inv:=1/lc_u mod p;

g,s,t:=eea_gcd(m_m,m_u,lc_m_inv,0,0,lc_u_inv,p);
Expand(s) mod p;
Expand(t) mod p;