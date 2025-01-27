# Black box for evaluating the rational function at a point
B:=proc(a,b,point_,p)
    local c:
    c:= a/b:
     print("c=",c):
    return [seq(Eval(c,x=i) mod p,i=point_[1..numelems(point_)])]:
end proc:
Eea_gcd:=proc(r0,r1,t0,t1,p)
    local r,t,q,i,f,g,qmax:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
     print("t0= ",t0):
     print("t1= ",t1):
    f:=r0:
    g:=t1:
    qmax:=1:
    i:=1:
    while r[i] <> 0 do
    #  1. find quotient and remainder
    q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p:
     print("r[",i,"]=",r[i]):
     print("q[",i,"]=",q[i]):
     print("t[",i,"]=",t[i]):
     print("Degree of q[",i,"]=",degree(q[i],x)):
     print("--------------------------------------"):
     print("r[",i+1,"]=",r[i+1]):
    if degree(q[i],x)>= qmax then 
        qmax:=degree(q[i],x):
        f:=r[i]:
        g:=t[i]:
    end if:
    #   reduction step to update t.
    t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
    #  Normalization step for r and t 
     if r[i+1]<> 0 then 
        rho:=lcoeff(r[i+1]):
        r[i+1]:=r[i+1]/rho mod p:
        t[i+1]:=t[i+1]/rho mod p:
    end if:
     print("f=",f):
     print("g=",g):
     print("______________________________________"):
    i:=i+1:
    end do:
    return f,g:# make g monic. Also add a check to see that g !=0
end proc:

p:=2^31-1:

num_deg:=7:
denom_deg:=5:
f0:=randpoly(x,degree=num_deg):
#  lc_f0:=lcoeff(f0):
#  f0:=f0/lc_f0 mod p:
f1:=randpoly(x,degree=denom_deg):
lc_f1:=lcoeff(f1):
f1:=f1/lc_f1 mod p:
# Black box find degree. 
d:=num_deg+denom_deg+1:
num_points:=d+1:
#  deg m = d > N+M
#  1. generate the points
points:=[seq(1..num_points)]:
y:=B(f0,f1,points,p):
m:=Expand(product(x-points[i],i=1..numelems(points))) mod p:
u:=Interp(points,y,x)mod p:# will be replaced by T(x,beta2*x-beta2*sig1+sig2) in seperation algorithm
# where sig=[sig1,sig2] is the evaluation vector and beta=[beta1,beta2] is a random vector.  
lc_u:=lcoeff(u):
lc_u_inv:=1/lc_u mod p:
m_u:=u/lc_u mod p:

#  2. generate the polynomial
f,g:=Eea_gcd(m,m_u,0,lc_u_inv,p):
lc_g:=lcoeff(g):
g:=g/lc_g mod p:
f:=f/lc_g mod p:
h:=f/g:
 print(h):
correct:=f0/f1 mod p:
print(simplify(correct)-simplify(h));
