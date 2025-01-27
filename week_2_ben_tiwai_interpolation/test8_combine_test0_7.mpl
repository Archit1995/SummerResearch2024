with(LinearAlgebra):
with(ArrayTools):
#  Black box
B:=proc(a,var,point_)
    local u,v;
    return [seq(eval(a,{seq(var[v]=point_[u][v],v=1..numelems(point_[u]))}),u=1..numelems(point_))];
end proc;
# Generating a prime for each variable 
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
# Generating a list of list powers of prime. 
generate_prime_powers:=proc(T,prime_points,num_var)
    local i,j;
    return [seq([seq(prime_points[j]^i,j = 1..num_var)], i = 0..2*T-1)];
end proc;  

vars:={x,y,z}:
f:=randpoly(vars,degree=5):
deg:=degree(f):
# Try T:=1..2^n until we find a T that works(O(log(n)) time complexity)
T:=deg+2:
num_var:=3;
prime_points:=generate_evaulation_primes(num_var);
prime_powers:=generate_prime_powers(T,prime_points,num_var);
y:=B(f,vars,prime_powers);