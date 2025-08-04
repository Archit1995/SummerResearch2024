MQRFR:=proc(r0,r1,t0,t1,p)
    print("In MQRFR"):
    local r,t,q,i,f,g,qmax,lcg:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
    f:=r0:
    g:=t1:
    qmax:=0:
    i:=1:
    while r[i] <> 0 do
    #  1. find quotient and remainder
    q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p:
    if degree(q[i],x)> qmax then 
        qmax:=degree(q[i],x):
        f:=r[i]:
        g:=t[i]:
    end if:
    t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
    i:=i+1:
    if qmax <=1 or gcd(f,g) <> 1 or g = 0 then 
        FAIL:
    end if:
    end do:
    lcg:=lcoeff(g):
    return f/lcg mod p,g/lcg mod p,qmax,lcg :
end proc:
