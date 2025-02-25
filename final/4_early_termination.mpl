early_termination_exponent_generator:=proc(T,prime_points,num_var,p)
    local i,j:
    print("In early_termination_exponent_generator");
    print("T=",T);
    print("prime_points=",prime_points);
    return modp([seq([seq(prime_points[j]^i,j = 1..num_var)], i = 0..2*T-1)],p):
end proc:  