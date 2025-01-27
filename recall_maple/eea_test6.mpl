eea:=proc(m,u)
    local r0,r1,t0,t1,f,g,qmax,i,r,t,q;
    if u=0 return (0,1) fi;
    r0:=m;
    r1:=u;
    t0:=0;
    t1:=1;
    f:=r1;
    g:=t1;
    qmax:=1;
    i:=1;
    while r[i]<>0 do
        q[i]:=quo(r[i-1],r[i]);
        if degree(q[i])>qmax then 
            qmax:=degree(q[i]) ;
            f:=r[i];
            g:=t[i];
        fi;
        r[i+1]:=r[i-1]-q[i]*r[i];
        t[i+1]:=t[i-1]-q[i]*t[i];
        i:=i+1;
    od;
    if qmax<=1 or gcd(f,g)<>1 then return FAIL fi;
    return (f/lcoeff(f),g/lcoeff(g));
end proc;