read "./TESTCASE_GENERATOR.mpl":
with(ArrayTools):
Construct_Blackbox:=proc(f,vars)
    local BB:
    BB:=proc(point_,p)
        local u,v;
        var:=vars:
        a:=f:
        return Eval(a,{seq(var[v]=point_[v],v=1..numelems(point_))}) mod p:
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
# 3. Generating a list of list powers of prime. [[2^0,3^0,5^0],[2^1,3^1,5^1],[2^2,3^2,5^2],..., [2^(2T-1),3^(2T-1),5^(2T-1)]]
generate_prime_powers:=proc(T,prime_points,num_var,p)
    local i,j:
    lprint("In generate_prime_powers");
    lprint("T=",T);
    lprint("prime_points=",prime_points);
    return modp([seq([seq(prime_points[j]^i,j = 1..num_var)], i = 0..2*T-1)],p):
end proc:  

# 4. get number of terms, lambda polynomial, and roots of lambda polynomial 
get_num_terms_lambda_roots:=proc(B,T,prime_points,num_var,p)# Needs correction find a way to keep roots in mem for t 
# and t+1 to check if the number of roots is equal to the number of terms. 
    lprint("In get_num_terms_lambda_roots");
    local t,flag,prime_powers,Y,Lambda,terms,R:
    t:=T:
    i:=1;
    while true do 
        lprint("i=",i):
        lprint("t=",t):
        prime_powers:=generate_prime_powers(t,prime_points,num_var,p):
        lprint("prime_powers=",prime_powers):
        Y:=[seq(B(prime_powers[p_p],p),p_p=1..numelems(prime_powers))]:
        lprint("Y=",Y):
        Lambda := BMEA(Y,p,Z):
        lprint("Lambda=",Lambda):
        terms[i]:=degree(Lambda,Z):
        lprint("terms=",terms[i]):
        R := Roots(Lambda) mod p:
        lprint("R=",R):
        lprint("num of roots of lambda=",nops(R)):
        # if i=1 then 
        if nops(R)<terms[i] or terms[i] <> terms [i-1]   then 
            t:=t*2:
        # else 
        # elif terms[i-1]=terms[i] and terms[i] = nops(R) then
        elif terms[i] = nops(R) and terms[i] = terms [i-1] then 
            lprint("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");
                return terms[i],Lambda,R,Y:
        end if: 
        # end if:
        i:=i+1:
        # if i = 5 then break end if:
    end do:
    lprint("terms=",terms[i]):
end proc:


BMEA := proc(v::list,p::posint,Z::name) # Might need correction because roots of lambda polynomial are not factors of 2,3,5 for all degrees
  local n,m,R0,R1,V0,V1,i:
  lprint("In BMEA");
  lprint("v=",v);    
  n := iquo( nops(v), 2 ):
  lprint("n=",n);
  m := 2*n-1:
  lprint("m=",m);
  R0 := Z^(2*n):
  lprint("R0=",R0);
  R1 := add( v[m+1-i]*Z^i, i=0..m ):
  lprint("R1=",R1);
# llprint("R0=",R0):
# llprint("R1=",R1):
  V0 := 0:
  V1 := 1:
  while n <= degree(R1,Z) do
     R0,R1 := R1,Rem(R0,R1,Z,'Q') mod p:
# llprint("R0=",R0):
# llprint("R1=",R1):
     V0,V1 := V1,Expand(V0-Q*V1) mod p:
# llprint("V0=",V0):
# llprint("V1=",V1):
  od:
  i := 1/lcoeff(V1,Z) mod p:
  return i*V1 mod p:
end:

generate_monomials:=proc(roots_,num_var,prime_points,vars)# needs correction- We are getting roots of lambda polynomial that are not 
    # factors of 2,3 and 5;
    local m,mm,i,j,counter,M_: 
    M_:=Vector(numelems(roots_),0):
    lprint("r=",numelems(roots_)):
    for i from 1 to  numelems(roots_) do # each root
        lprint("i=",i):
        mm:=roots_[i]:
        m:=1:
        lprint("roots_[i]=",roots_[i]):
        for j from 1 to numelems(prime_points) do #  each prime
            counter:=0:
            lprint("j=",j):
            while mm mod prime_points[j] = 0 do #repeated division
                lprint("prime_points[j]=",prime_points[j]):
                mm:=iquo(mm,prime_points[j]):
                lprint("mm=",mm):
                counter:=counter+1:
                lprint("counter=",counter):
                lprint("================================================"):
            end do:
            m:=m*vars[j]^counter:# each monomial
            lprint("m=",m):
            lprint("-----------------------------------------------------------------"):
        end do:
        lprint("m=",m):
        lprint("i=",i):
        M_[i]:=m:
        lprint("M[i]=",M_[i]):
        lprint("______________________________________"):
    end do:
    lprint([seq(ifactor(roots_[i]),i=1..numelems(roots_))]);# We are getting roots of lambda polynomial that are not 
    # factors of 2,3 and 5;
    return convert(M_,list):
end proc:


Zippel_Vandermonde_solver:=proc(Y_::list,terms::integer,roots_::list,lambda_,p::integer)# Correct so far. 
    local M,fin_coeff,q,q_lambda_inv,V_inv_b,i,j,y:
    M:=lambda_ mod p:
    lprint("In Zippel_Vandermonde_solver"):
    y:=Y_ mod p:
    lprint("y=",y):
    lprint("terms=",terms):
    lprint("roots_=",roots_):
    lprint("lambda_=",lambda_):
    lprint("M=",M):
    lprint("roots=",roots_):
    fin_coeff:=Vector(terms,0):
    for i from 1 to terms do
        q:=quo(M,Z-roots_[i],Z) mod p:
        lprint("q=",q):
        q_lambda_inv:= 1/ Eval(q,Z=roots_[i]) mod p:
        lprint("q_lambda_inv=",q_lambda_inv):
        qq:=q*q_lambda_inv mod p;
        lprint("q*q_lambda_inv mod p=",qq):
        V_inv_b:=0:
        for j from 1 to terms do
            lprint("j=",j):
            lprint("coeff(q,Z,j-1)=",coeff(q,Z,j-1)):
            lprint("y[j]=",y[j]):
            
            lprint("coeff(q,Z,j-1)*y[j]=",coeff(q,Z,j-1)*y[j] mod p):
            # lprint("__________________________________");
            # lprint("coeff(qq,Z,j-1)=",coeff(qq,Z,j-1)):
            # lprint("coeff(qq,Z,j-1)*y[j]=",coeff(qq,Z,j-1)*y[j] mod p):
            # lprint("__________________________________");
            V_inv_b:=V_inv_b+coeff(q,Z,j-1)*y[j] mod p:
            lprint("V_inv_b after adding term for j=",j," is ",V_inv_b):
        end do:
        lprint("__________________________________");
        lprint("V_inv_b=",V_inv_b):
        lprint("q_lambda_inv=",q_lambda_inv):
        lprint("__________________________________");
        fin_coeff[i]:=V_inv_b*q_lambda_inv mod p:
        lprint("fin_coeff[",i,"]=",fin_coeff[i]):
        lprint("--------------------------------------");
    end do:
    lprint("fin_coeff=",fin_coeff):
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

# test_case:="bspline_small_sys_low_deg2":

test_case:=2:

vars,p,T,ff,gg:=data_generator(test_case):
T:=4:
p:=17:
num_var:=nops(vars):
# num_var:=30:
# num_terms:=32:
# den_terms:=44:
# num_var:=21:
# num_terms:=1033:
# den_terms:=11:
# vars,p,T,ff,gg:=data_generator(test_case,num_var,num_terms,den_terms):
prime_points:=generate_evaulation_primes(num_var):
f:=ff:
# f:=gg:
lprint("f= ",f):
lprint("vars= ",vars):
B:=Construct_Blackbox(f,vars);

terms,Lambda,R,Y:=get_num_terms_lambda_roots(B,T,prime_points,num_var,p):

lprint("R=",R):
Roots_ := [ seq(r[1], r in R ) ]:
lprint("Roots_=",Roots_):
Monomials:=generate_monomials(Roots_,num_var,prime_points,vars):
lprint("Monomials=",Monomials):
coeff_:= Zippel_Vandermonde_solver(Y,terms,Roots_,Lambda,p):
lprint("coeff_=",coeff_):    
# a:=-62*x^2*z^3+97*x*y^3*z-73*y*z^4-56*x*y*z^2 +87*x*y mod p:
# a:=2*x*y+3*z+1 mod p:
f;
f1:=construct_final_polynomial(coeff_,Monomials);
f1-f mod p;
# g:=x+y:
# Bg:=Construct_Blackbox(g,vars);

# terms_g,Lambda_g,R_g,Y_g:=get_num_terms_lambda_roots(Bg,T,prime_points,num_var,p):

# Roots_g := [ seq(r[1], r in R_g ) ]:
# Monomials_g:=generate_monomials(Roots_g,num_var,prime_points,vars):
# coeff_g:= Zippel_Vandermonde_solver(Y_g,terms_g,Roots_g,Lambda_g,p):
# g1:=construct_final_polynomial(coeff_g,Monomials_g);
# g1;
# g1-g;