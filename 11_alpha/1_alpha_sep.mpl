# 1. Have a black box for f/g. Change black box implementation- it only takes in the point and the prime.
Construct_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        local u,v,var,num,den,a;
        var:=vars:
        num:=f:
        denom_:=g:
        a:=num/denom_:
        # print("a: ",a);
        return Eval(a,{seq(var[i]=point_[i],i=1..numelems(point_))}) mod p: 
    end proc:
    return BB:
end proc:
 generate_random_vector:=proc(n,p)
    r:=rand(p);
    return [seq(r(),i=1..n)];
end proc:

# MQRFR(m,u,0,1,p)
MQRFR:=proc(r0,r1,t0,t1,p)
    local r,t,q,i,f,g,qmax:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
    #  print("t0= ",t0):
    #  print("t1= ",t1):
    f:=r0:
    g:=t1:
    qmax:=0:
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
    #  if r[i+1]<> 0 then 
    #     rho:=lcoeff(r[i+1]):
    #     r[i+1]:=r[i+1]/rho mod p:
    #     t[i+1]:=t[i+1]/rho mod p:
    # end if:
     print("f=",f):
     print("g=",g):
     print("______________________________________"):
    i:=i+1:
    if qmax <=1 or gcd(f,g) <> 1 then 
        FAIL:
    end if:
    end do:
    lcg:=lcoeff(g):
    return f/lcg mod p,g/lcg mod p,qmax:# make g monic. Also add a check to see that g !=0
end proc:

vars:={x,y}:
num_var:=numelems(vars):
p:=2^31-1:
num_deg:=3:
denom_deg:=2:
ff:=randpoly(vars,sparse,degree=num_deg) mod p;
gg:=randpoly(vars,sparse,degree=denom_deg) mod p;
print("lcoeff(gg)= ",lcoeff(gg) mod p);
gg:=gg/lcoeff(gg) mod p;
print("ff ",ff);
print("gg ",gg);
B:=Construct_Blackbox(ff,gg,vars);
print(B);

num_points:=1:
correct_degree:=true:
while(correct_degree) do
    print("num_points: ",num_points);
    r:=rand(p):
    alpha:=[seq(r(),i=1..num_points)];
    print("alpha: ",alpha);
    # sigma[numpoints][num_var]
    points:=[seq(generate_random_vector(numelems(vars),p),i=1..num_points)];
    print("points: ",points);
    # beta[numpoints][num_var-1]
    shift_:=[seq(generate_random_vector(numelems(vars)-1,p),i=1..num_points)];# beta
    print("shift_: ",shift_);                                                                                     
    # Projection phi(x,beta[numpoints][numvar-1] *alpha[numpoints]-beta[numpoints][numvar-1]*sigma[numpoints][1]+
    # sigma[numpoints][numvar])
    phi_:=[seq([seq(0,nv=1..num_var)],np=1..num_points)];
    # print("phi_: ",phi_);
    # phi_[np][1]=alpha[np];
    for np from 1 to num_points do 
        phi_[np][1]:=alpha[np];
        for nv from 2 to num_var do 
            phi_[np][nv]:=shift_[np][nv-1]*alpha[np]-shift_[np][nv-1]*points[np][1]+points[np][nv] mod p;
        end do;
    end do;
    print("phi_: ",phi_);
#   phi_[np][nv]=beta[np][nv-1]*alpha[np]-beta[np][nv-1]*sigma[np][1]+sigma[np][nv];
    # point_:=[seq(shift__[2]*x-shift__[2]*point__[1]+point__[2];)]
    Y:=[seq(B(phi_[i],p),i=1..num_points)];
    # print("Y: ",Y);
    m:=Expand(product(x-alpha[i],i=1..numelems(points))) mod p;
    print("m: ",m);
    # points_for_uni_interp:=[seq(points[i][1],i=1..num_points)];
    # print("points_for_uni_interp: ",points_for_uni_interp);
    u:=Interp(alpha,Y,x)mod p:
     print("u: ",u);
    checking:=[seq(Eval(u,x=alpha[i]) mod p,i=1..num_points)];
    print("Y: ",Y);
    print("checking: ",checking);
   
    # lc_u:=lcoeff(u):
    # lc_u_inv:=1/lc_u mod p:
    # m_u:=u/lc_u mod p:
    f,g,dq:=MQRFR(m,u,0,1,p):
    # if num_points = 10 then break; end if;
    if dq > 1 then 
        break:
    else 
        num_points:=num_points*2:
    end if;
    if num_points > 16 then break; end if;
    # end if:
    print("====================================================="):
end do:

# lc_g:=lcoeff(g):
# g:=g/lc_g mod p:
# f:=f/lc_g mod p:
# h:=f/g:
#  print(h):
# correct:=ff/gg mod p:
# print("correct = ",correct);
# print(simplify(correct)-simplify(h));



# sigma_:=[1,2];
# beta_:=[3,5];# Randomize betas so that we can use sigma as primes for monomial evaluation at primes. 
# beta_[2]*x-beta_[2]*sigma_[1]+sigma_[2];

shift__[2]*x-shift__[2]*point__[1]+point__[2];
# temp:=B([x,beta_[2]*x-beta_[2]*sigma_[1]+sigma_[2]],p);
# T:=modp(simplify(temp[1]),p)/modp(simplify(temp[2]),p); 
# print(T);
# In mqrfr we need m to be 
#  m:=Expand(product(x-points[i],i=1..numelems(points))) mod p:
# But points are vectors here, we need one dimensional points. 

# T:=num/denom 
# Deg(num)=degree(num,x)+degree(num,y);
# Deg(denom)=degree(denom,x)+degree(denom,y);
# Make the denominator monic. And dividide the numerator by the leading coefficient of the denominator.