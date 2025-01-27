test_sum:=proc(a,b)
    return a+b;
end proc;
test_mod:=proc(a,b)
    local q,r;
    q:= quo(a,b,x,'r');
    return r;
end proc;
f1:=randpoly({x,y},degree=2);
f2:=randpoly({x,y},degree=3);
test_sum(f1,f2);
test_mod(f2,f1);
