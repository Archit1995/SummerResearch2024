# 1. Have a black box for f/g. Change black box implementation- it only takes in the point and the prime.
Construct_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        local u,v,var,num,den,a;
        var:=vars:
        num:=f:
        denom_:=g:
        a:=num/denom_:
        # print("a: ",a);
        return Eval(a,{seq(var[i]=point_[i],i=1..numelems(point_))}) mod p: 
    end proc:
    return BB:
end proc:
 

vars:={x,y}:
# prime_points:=generate_evaulation_primes(num_var):
# prime_points:
# p:=2^31-1:
p:=19;
f:=randpoly(vars,sparse,degree=7) mod p:
g:=randpoly(vars,sparse,degree=5) mod p:
B:=Construct_Blackbox(f,g,vars);
print(B);
B([1,2],p);
# sigma_:=[1,2];
# beta_:=[3,5];# Randomize betas so that we can use sigma as primes for monomial evaluation at primes. 
# beta_[2]*x-beta_[2]*sigma_[1]+sigma_[2];
# temp:=B([x,beta_[2]*x-beta_[2]*sigma_[1]+sigma_[2]],p);
# T:=modp(simplify(temp[1]),p)/modp(simplify(temp[2]),p); 
# print(T);
# In mqrfr we need m to be 
#  m:=Expand(product(x-points[i],i=1..numelems(points))) mod p:
# But points are vectors here, we need one dimensional points. 

# T:=num/denom 
# Deg(num)=degree(num,x)+degree(num,y);
# Deg(denom)=degree(denom,x)+degree(denom,y);
# Make the denominator monic. And dividide the numerator by the leading coefficient of the denominator.