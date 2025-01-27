with(LinearAlgebra):
with(ArrayTools):
# 1. Black box for some polynomial f in Z[x_1,x_2,....x_n] of some degree m
B:=proc(point_,p)
    local u,v,a,var:
    var:={x,y,z}:
    a:=-62*x^2*z^3+97*x*y^3*z-73*y*z^4-56*x*y*z^2 +87*x*y mod p:
    return modp([seq(eval(a,{seq(var[v]=point_[u][v],v=1..numelems(point_[u]))}),u=1..numelems(point_))],p):
end proc:
# 2. Generating a prime for each variable 
generate_evaulation_primes:=proc(n)
    local p,m,i:
    m:=1:
    p:=Vector(n,0):
    for i from 1 to n do 
        p[i]:=nextprime(m):
        m:=p[i]:
    end do:
return convert(p,list):
end proc:   
# 3. Generating a list of list powers of prime. 
generate_prime_powers:=proc(T,prime_points,num_var)
    local i,j:
    return [seq([seq(prime_points[j]^i,j = 1..num_var)], i = 0..2*T-1)]:
end proc:  

# 4. Getting the number of terms in the polynomial
get_num_terms:=proc(v,T)
    local H,i:
    H:=Matrix([seq(v[i..i+(T-1)],i=1..T)]):
    return H,Rank(H):
end proc:
# 5. Getting the roots of the lambda polynomial
get_rootsOf_lambda_polynomial:=proc(M,v,terms)
    local H,b,X,num_row,Lambda,R,i,r:
    H:=M[1..terms,1..terms]:
    b:=-Vector(v[terms+1..terms+terms]):
    X:=LinearSolve(H,b):
    num_row:=Size(X)[1]:
    Lambda:=Z^num_row:
    for i from 1 to num_row do
        Lambda:=Lambda+X[i]*Z^(i-1):
    end do:
    R:=roots(Lambda):
    return [seq(r[1],r in R)]:
end proc:

# 6.Generating monomials from the roots of the lambda polynomial
generate_monomials:=proc(roots_,num_var,prime_points,vars)
     local ff,l,l2,i,prime_var_map,monomials,j:
     prime_var_map:= table([seq(prime_points[i]=vars[i],i=1..num_var)]):
     monomials:=Vector(numelems(roots_),0):
     for j from 1 to numelems(roots_) do 
        #   print(roots_[j]):
          ff:=ifactor(roots_[j]):
        #   print(ff):
          l:=nops(ff):
          for i from 1 to l do
               l2:=nops(op(i,ff)):
               if l2=1 then 
                    ff:=subs(op(i,ff)=prime_var_map[op(1,op(i,ff))],ff):
               else 
                    ff:=subs(op(1,op(i,ff))=prime_var_map[op(1,(op(1,op(i,ff))))],ff):
               fi:
          end do:
          monomials[j]:=ff:
     end do:
     return convert(monomials,list):
end proc:

# Step 2 of BT interpolation
# 7. Constructing the Vandermonde matrix - Replace with Zippel Vandermonde
Construct_Vandermonde:=proc(terms,Roots_)
    local i,j:
    return Matrix([seq([seq(Roots_[j]^i,j = 1..numelems(Roots_))], i = 0..terms-1)]):    
end proc:

# 8. Getting the coefficients of the polynomial
get_coefficients:=proc(terms,Roots_,v)
    local Van,b:
    b:=<v[1..terms]>:
    Van:=Construct_Vandermonde(terms,Roots_):
    return LinearSolve(Van,b):
end proc:
# 9. Constructing the final polynomial
construct_final_polynomial:=proc(coeff_,Monomials)
    local i,f,n:
    f:=0:
    for i from 1 to numelems(coeff_) do
        f:=f+coeff_[i]*Monomials[i]:
    end do:
    return f:
end proc:


get_num_terms:=proc(prime_points,num_var,p)
    local T,H,i,v,Y,prime_powers,term:
    i:=2:
    term[i]:=-2:
    term[i-1]:=-1:
    T:=2:
    while term[i-1]<>term[i]do 
        prime_powers:=generate_prime_powers(T,prime_points,num_var):
        v:=B(prime_powers,p):
         H:=Matrix([seq(v[i..i+(T-1)],i=1..T)]):
        i:=i+1:
        T:=T*2:
        term[i]:=Rank(H):
    end do:
    return term[i],H,v:
end proc:

prime_points:

# Tester 
p:=7;
num_var:=3:
vars:={x,y,z}:
prime_points:=generate_evaulation_primes(num_var):

terms,Y,y_:=get_num_terms(prime_points,num_var,p):
Roots_:=get_rootsOf_lambda_polynomial(Y,y_,terms):
Monomials:=generate_monomials(Roots_,num_var,prime_points,vars):
coeff_:=get_coefficients(terms,Roots_,y_):
f1:=construct_final_polynomial(coeff_,Monomials);
