
generate_evaulation_primes:=proc(n)
    local p,m,i;
    m:=1;
    p:=Vector(n,0);
    for i from 1 to n do 
        p[i]:=nextprime(m);
        m:=p[i];
    end do;
return convert(p,list);
end proc;   
    