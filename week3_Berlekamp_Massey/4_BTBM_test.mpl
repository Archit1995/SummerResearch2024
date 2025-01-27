with(ArrayTools):
Construct_Blackbox:=proc(f,vars)
    local BB:
    BB:=proc(point_,p)
        local u,v;
        var:=vars:
        a:=f:
        return modp([seq(eval(a,{seq(var[v]=point_[u][v],v=1..numelems(point_[u]))}),u=1..numelems(point_))],p):
    end proc:
    return BB:
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
generate_prime_powers:=proc(T,prime_points,num_var,p)
    local i,j:
    return modp([seq([seq(prime_points[j]^i,j = 1..num_var)], i = 0..2*T-1)],p):
end proc:  

# 4. get number of terms, lambda polynomial, and roots of lambda polynomial 
get_num_terms_lambda_roots:=proc(B,T,prime_points,num_var,p)
    local t,flag,prime_powers,Y,Lambda,terms,R:
    t:=T:
    flag:=true:
    while flag do 
        # print("t=",t):
        prime_powers:=generate_prime_powers(t,prime_points,num_var,p):
        # print("prime_powers=",prime_powers):
        Y:=B(prime_powers,p):
        # print("Y=",Y):
        Lambda := BMEA(Y,p,Z):
        # print("Lambda=",Lambda):
        terms:=degree(Lambda,Z):
        # print("terms=",terms):
        R := Roots(Lambda) mod p:
        # print("R=",R):
        # print("num of roots of lambda=",nops(R)):
        if nops(R)<terms then 
            t:=t*2:
        end if:
        if nops(R)=terms then 
            flag:=false:
        end if:
    end do:
    return terms,Lambda,R,Y:
end proc:


BMEA := proc(v::list,p::posint,Z::name) 
  local n,m,R0,R1,V0,V1,i:
  n := iquo( nops(v), 2 ):
  m := 2*n-1:
  R0 := Z^(2*n):
  R1 := add( v[m+1-i]*Z^i, i=0..m ):
# lprint("R0=",R0):
# lprint("R1=",R1):
  V0 := 0:
  V1 := 1:
  while n <= degree(R1,Z) do
     R0,R1 := R1,Rem(R0,R1,Z,'Q') mod p:
# lprint("R0=",R0):
# lprint("R1=",R1):
     V0,V1 := V1,Expand(V0-Q*V1) mod p:
# lprint("V0=",V0):
# lprint("V1=",V1):
  od:
  i := 1/lcoeff(V1,Z) mod p:
  i*V1 mod p:
end:

generate_monomials:=proc(roots_,num_var,prime_points,vars)
    local m,mm,i,j,counter,M_: 
    M_:=Vector(numelems(roots_),0):
    # print("r=",numelems(roots_)):
    for i from 1 to  numelems(roots_) do # each root
        # print("i=",i):
        mm:=roots_[i]:
        m:=1:
        # print("roots_[i]=",roots_[i]):
        for j from 1 to numelems(prime_points) do #  each prime
            counter:=0:
            # print("j=",j):
            while mm mod prime_points[j] = 0 do #repeated division
                # print("prime_points[j]=",prime_points[j]):
                mm:=iquo(mm,prime_points[j]):
                # print("mm=",mm):
                counter:=counter+1:
                # print("counter=",counter):
                # print("================================================"):
            end do:
            m:=m*vars[j]^counter:# each monomial
            # print("m=",m):
            # print("-----------------------------------------------------------------"):
        end do:
        # print("m=",m):
        # print("i=",i):
        M_[i]:=m:
        # print("M[i]=",M_[i]):
        # print("______________________________________"):
    end do:
    return convert(M_,list):
end proc:


Zippel_Vandermonde_solver:=proc(y::list,terms::integer,roots_::list,lambda_,p::integer)
    local M,fin_coeff,q,q_lambda_inv,V_inv_b,i,j:
    M:=lambda_ mod p:
    # print("M=",M):
    # print("roots=",roots_):
    fin_coeff:=Vector(terms,0):
    for i from 1 to terms do
        q:=quo(M,Z-roots_[i],Z):
        # print("q=",q):
        q_lambda_inv:= 1/ Eval(q,Z=roots_[i]) mod p:
        # print("q_lambda_inv=",q_lambda_inv):
        V_inv_b:=0:
        for j from 1 to terms do
            V_inv_b:=V_inv_b+coeff(q,Z,j-1)*y[j] mod p:
        end do:
        fin_coeff[i]:=V_inv_b*q_lambda_inv mod p:
    end do:
    # print("fin_coeff=",fin_coeff):
    return convert(fin_coeff,list):
end proc:

construct_final_polynomial:=proc(coeff_,Monomials)
    local i,f,n:
    f:=0:
    for i from 1 to numelems(coeff_) do
        f:=f+coeff_[i]*Monomials[i]:
    end do:
    return f:
end proc:

test:=proc(degrees)
    local num_var,vars,prime_points,p,f,B,T,terms,Lambda,R,Y,Roots_,Monomials,coeff_,f1,i,flag,a,counter,results,d,r:
    num_var:=3:
    vars:={x,y,z}:
    prime_points:=generate_evaulation_primes(num_var):
    prime_points:
    p:=2^31-1:
    results:=Vector(numelems(degrees),0):
    for d from 1 to  numelems(degrees) do
        print("degree=",degrees[d]); 
        f:=randpoly(vars,degree=degrees[d]) mod p:
        B:=Construct_Blackbox(f,vars);
        print(B);
        T:=2:
        terms,Lambda,R,Y:=get_num_terms_lambda_roots(B,T,prime_points,num_var,p):
        Roots_ := [ seq(r[1], r in R ) ]:
        Monomials:=generate_monomials(Roots_,num_var,prime_points,vars):
        coeff_:= Zippel_Vandermonde_solver(Y,terms,Roots_,Lambda,p):
        f1:=construct_final_polynomial(coeff_,Monomials);
        results[d]:=f1-f;
        print("Results=",results[d]);
        print("______________________________________________");
    end do:
    return results;
end proc:
init_:=5;
degrees:=[seq(i,i=init_..20)];
print("numelems(degrees",numelems(degrees));
Result:=test(degrees):
c:=0;
for rr from 1 to numelems(Result) do 
    if Result[rr]<>0 then 
        c:=c+1;
        print("degree = ",init_+rr-1);
        print(Result[rr]);
    end if:
end do:
print("number of failed tests=",c);
print("number of tests=",numelems(Result));
print(evalf(c/numelems(Result))*100);