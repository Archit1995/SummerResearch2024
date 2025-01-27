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
    eea_gcd(r1,r2,s1,t1,s2,t2);
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
f1:=randpoly(x,degree=num_deg);

f2:=randpoly(x,degree=denom_deg);
# y:=B(f1,f2,points);
# u:=interp(y,points,x);

g,s,t:=eea_gcd(f1,f2,1,0,0,1);
g-expand(s*f1+t*f2);