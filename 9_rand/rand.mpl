p:=2^32-1;
generate_random_vector:=proc(n,p)
    r:=rand(p);
    return [seq(r(),i=1..n)];
end proc:
# rr:=generate_random_vector(2,p);
rr:=[seq(generate_random_vector(2,p),i=1..10)];
print("rr= ",rr);