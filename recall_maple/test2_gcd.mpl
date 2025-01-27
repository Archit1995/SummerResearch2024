test_gcd:=proc(r0,r1)
    q:=quo(r0,r1,x,'r');
    print("r = ",r);
    if r = 0 then 
        return r1/lcoeff(r1);
    else
        test_gcd(r1,r);
    end if;
end proc;
f1:=randpoly({x,y},degree=3);
f2:=expand(f1*(x-1));
test_gcd(f2,f1);
ff:=expand((x-1)*(x-2)*(x-3)*(x-5)*(x-6)*(x-7));
ff2:=expand((x-1)*(x-11));
simplify(test_gcd(ff,ff2));