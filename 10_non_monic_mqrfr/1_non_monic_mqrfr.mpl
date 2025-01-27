# Black box for evaluating the rational function at a point
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
    #  print("r[",i,"]=",r[i]):
    #  print("q[",i,"]=",q[i]):
    #  print("t[",i,"]=",t[i]):
     print("Degree of q[",i,"]=",degree(q[i],x)):
     print("--------------------------------------"):
    #  print("r[",i+1,"]=",r[i+1]):
    if degree(q[i],x)>= qmax then 
        qmax:=degree(q[i],x):
        f:=r[i]:
        g:=t[i]:
        print("r[",i-1,"]=",r[i-1]):
        print("q[",i,"]=",q[i]):
        print("f=",f):
        print("g=",g):
    end if:
    #   reduction step to update t.
    t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
    #  Normalization step for r and t 
    #  if r[i+1]<> 0 then 
    #     rho:=lcoeff(r[i+1]):
    #     r[i+1]:=r[i+1]/rho mod p:
    #     t[i+1]:=t[i+1]/rho mod p:
    # end if:
    #  print("f=",f):
    #  print("g=",g):
    #  print("______________________________________"):
    i:=i+1:
    if qmax <=1 or gcd(f,g) <> 1 then 
        FAIL:
    end if:
    end do:
    lcg:=lcoeff(g):
    return f/lcg mod p,g/lcg mod p,qmax:# make g monic. Also add a check to see that g !=0
end proc:

# Black box find degree. 
Early_termination_MQRFR:=proc(B,p)
    local r,u,f,g,dq,num_points,correct_degree,points,Y,m:
    r:=rand(p):
    num_points:=1:
    correct_degree:=true:
    while(correct_degree) do 
        #  deg m = d > N+M
        #  1. generate the points
        # point_:=[seq(1..num_points)]:
        print("num_points= ",num_points):
        points:=[seq(r(),i=1..num_points)]:
        # points:=[seq(i,i=1..num_points)]:
        # print("points= ",points):
        # break:
        Y:=[seq(B(points[i],p),i=1..num_points)]:
        # y:=[seq(B(i,p),i=points)]:
        # print("y= ",y):
        m:=Expand(product(x-points[i],i=1..numelems(points))) mod p:
        print("m= ",m):
        u:=Interp(points,Y,x)mod p:# will be replaced by T(x,beta2*x-beta2*sig1+sig2) in seperation algorithm
        # where sig=[sig1,sig2] is the evaluation vector and beta=[beta1,beta2] is a random vector.  
        # print("u= ",u):
        if u = 0 then 
            f:=0;
            g:=1;
            break:
        else 
            f,g,dq:=MQRFR(m,u,0,1,p):
        end if:

        if dq > 1 then 
            break:
        else 
            num_points:=num_points*2:
        end if:
        # print("====================================================="):
    end do:
    return f,g
end proc:



p:=2^31-1:
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
print(B);
f,g:=Early_termination_MQRFR(B,p):
print("f= ",f):
print("g= ",g):
h:=f/g:
 print(h):
correct:=f0/f1 mod p:
print(simplify(correct)-simplify(h));
