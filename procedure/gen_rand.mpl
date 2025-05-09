generate_random_vectors:=proc(n,p)
    local i,r:
    r:=rand(p):
    return [seq(r(),i=1..n)];
end proc: