exp_:=proc(v,l,p)
    return([seq(v[i]^l mod p,i=1..nops(v))]);
end proc:
