B:=proc(a,b,point_,p)
    local c;
    c:= a/b;
    print("c=",c);
    return [seq(Eval(c,x=i) mod p,i=point_[1..numelems(point_)])];
end proc;
MEea_gcd:=proc(r0,r1,t0,t1,p,N)
    print("r0=",r0);
    print("t0=",t0);
   local r2,q,s2,t2;
    if r1 = 0 or degree(r1,x)<= N then return r1, t1 end if;
    # 1. find quotient and remainder
    q:= Quo(r0,r1,x,'r2') mod p;
    print("q=",q);
    print("r2=",r2);
    #  reduction step to update s and t.
    # s2:=s0-q*s1;
    t2:=Expand(t0-q*t1)mod p:
    # t2:=t2/lcoeff(t2) mod p ;
    print("t2=",t2);
    # Normalization
    if r2<> 0 then 
        rho:=lcoeff(r2);
        r2:=r2/rho mod p;
        t2:=t2/rho mod p;
    end if;
    print("-----------------");
    MEea_gcd(r1,r2,t1,t2,p,N);
end proc;
p:=19;

f0:=3*x^2+4;
f1:=x+2;

num_deg:=2;
denom_deg:=1;
d:=num_deg+denom_deg;
num_points:=d+1;
# deg m = d > N+M
# 1. generate the points
points:=[15,13,6,2];

y:=B(f0,f1,points,p);
# y:=[13, 10, 14, 4];

u:=Interp(points,y,x)mod p;
lc_u:=lcoeff(u);
lc_u_inv:=1/lc_u mod p;
m_u:=u/lc_u mod p;

m:=Expand(product(x-points[i],i=1..numelems(points))) mod p;
lc_m:=lcoeff(m);
m_m:=m/lc_m mod p;
# 2. generate the polynomial

f,g:=MEea_gcd(m_m,m_u,0,lc_u_inv,p,num_deg);
lc_g:=lcoeff(g);
g:=g/lc_g mod p;
f:=f/lc_g mod p;
# g-expand(s*f1+t*f2);
print(f/g);