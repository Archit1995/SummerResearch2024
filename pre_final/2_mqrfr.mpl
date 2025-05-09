MQRFR:=proc(r0,r1,t0,t1,p)
    #print("In MQRFR"):
    local r,t,q,i,f,g,qmax,lcg:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
     #print("r0= ",r0):
     #print("t1= ",t1):
    f:=r0:
    g:=t1:
    qmax:=0:
    i:=1:
    while r[i] <> 0 do
    #  1. find quotient and remainder
    #print("i= ",i):
    q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p:

     #print("r[",i-1,"]=",r[i-1]):
     #print("r[",i,"]=",r[i]):
     #print("r[",i+1,"]=",r[i+1]):
     #print("q[",i,"]=",q[i]):
     #print("Degree of q[",i,"]=",degree(q[i],x)):
     #print("t[",i,"]=",t[i]):
     
     
    #  #print("--------------------------------------"):
     
    if degree(q[i],x)> qmax then 
        qmax:=degree(q[i],x):
        f:=r[i]:
        g:=t[i]:
    end if:
    #   reduction step to update t.
    #print("lead coeff of g=",lcoeff(g)):
    t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
    #print("t[",i+1,"]=",t[i+1]):
    #  Normalization step for r and t 
    #  if r[i+1]<> 0 then 
    #     rho:=lcoeff(r[i+1]):
    #     r[i+1]:=r[i+1]/rho mod p:
    #     t[i+1]:=t[i+1]/rho mod p:
    # end if:
    #  #print("f=",f):
    #  #print("g=",g):
    #print("______________________________________"):
    i:=i+1:
    if qmax <=1 or gcd(f,g) <> 1 or g = 0 then 
        FAIL:
    end if:
    # break:
    end do:
    lcg:=lcoeff(g):
    # #print("lcg= ",lcg):
    return f/lcg mod p,g/lcg mod p,qmax,lcg :# make g monic. Also add a check to see that g !=0
end proc:

# With normailzation step sum of lcoeff(g)=1/lcg mod p ??