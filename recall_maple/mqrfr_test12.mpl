B:=proc(a,b,point_,p)
    local c;
    c:= a/b;
    print("c=",c);
    return [seq(Eval(c,x=i) mod p,i=point_[1..numelems(point_)])];
end proc;
Eea_gcd:=proc(r0,r1,t0,t1,p)
    local r,t,q,i,f,g,qmax;
    r[0]:=r0;
    r[1]:=r1;
    t[0]:=t0;
    t[1]:=t1;
    print("t0= ",t0);
    print("t1= ",t1);
    f:=r0;
    g:=t1;
    qmax:=1;
    i:=1;
    while r[i] <> 0 do
    # 1. find quotient and remainder
    q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p;
    print("r[",i,"]=",r[i]);
    print("q[",i,"]=",q[i]);
    print("t[",i,"]=",t[i]);
    print("Degree of q[",i,"]=",degree(q[i],x));
    print("--------------------------------------");

    # print("r[",i+1,"]=",r[i+1]);
    if degree(q[i],x)>= qmax then 
        qmax:=degree(q[i],x);
        f:=r[i];
        g:=t[i];
    end if;
    #  reduction step to update t.
    t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
    # Normalization step for r and t 
     if r[i+1]<> 0 then 
        rho:=lcoeff(r[i+1]);
        r[i+1]:=r[i+1]/rho mod p;
        t[i+1]:=t[i+1]/rho mod p;
    end if;
    print("f=",f);
    print("g=",g);
    print("______________________________________");
    i:=i+1;
    end do;
    # if dq <=1 or gcd(r2,t2)<>1 then return FAIL end if;
    # Eea_gcd(f,r2,g,t2,p,dq);
    return f,g;
end proc;

p:=19;

f0:=3*x^2+4;
f1:=x+2;

num_deg:=2;
denom_deg:=1;
d:=num_deg+denom_deg+1;
num_points:=d+1;
# deg m = d > N+M
# 1. generate the points
points:=[0,1,7,2,13];

y:=B(f0,f1,points,p);
# y:=[13, 10, 14, 4];

u:=Interp(points,y,x)mod p;
lc_u:=lcoeff(u);
lc_u_inv:=1/lc_u mod p;
m_u:=u/lc_u mod p;

m:=Expand(product(x-points[i],i=1..numelems(points))) mod p;

# 2. generate the polynomial

f,g:=Eea_gcd(m,m_u,0,lc_u_inv,p);
lc_g:=lcoeff(g);
g:=g/lc_g mod p;
f:=f/lc_g mod p;
# g-expand(s*f1+t*f2);
print(f/g);