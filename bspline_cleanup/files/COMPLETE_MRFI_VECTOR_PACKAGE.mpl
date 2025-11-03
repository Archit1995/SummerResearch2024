########################################################
# COMPLETE MULTIVARIATE RATIONAL FUNCTION INTERPOLATION
# FOR VECTOR-VALUED BLACK BOXES (LINEAR SYSTEMS)
########################################################
# This file contains ALL necessary procedures in one place
# for easier testing and debugging

print("Loading Complete MRFI Vector Package..."):

########################################################
# 0. CONSTRUCT BLACK BOX - For Linear Systems
########################################################

Constuct_Sys_Blackbox:=proc(System,Vars,params) 
    local Lin_BB:
    Lin_BB:=proc(point_,p)option remember:
        local unordered_soln,Soln,soln,var,i,v:
        print("point_",point_):
        var:=params:
        Soln:=get_eqn(Sys,Vars):
        unordered_soln:=convert(Soln,list):
        soln:=reording(unordered_soln,nops(Sys)):
        return [seq(eval(op(2,soln[i]),{seq(var[v]=point_[v],v=1..numelems(point_))}),i=1..nops(Sys))] mod p:
    end proc:
    return Lin_BB:
end proc:

# Helper for solving the system
get_eqn:=proc(Sys,vars)option remember:
    print("in get_eqn"):
    return solve(Sys,Vars):
end proc:

########################################################
# 1. REORDERING HELPER FUNCTIONS
########################################################

reording:=proc(unordered_soln,num_eqn)
    local temp_var,component1,ordered_soln,i;
    ordered_soln:=[seq(0,i=1..num_eqn)];
    for i from 1 to num_eqn do 
        component1:=get_component(op(1,unordered_soln[i])):
        ordered_soln[component1]:=unordered_soln[i];
    end do;
    return ordered_soln;
end proc:

get_component:=proc(expression)
    local temp_var;
    temp_var:=convert(expression,string);
    return parse(temp_var[2..length(temp_var)]);
end proc:

########################################################
# 2. PROJECTION PHI
########################################################

projection_image_phi:=proc(num_var,alpha,beta_,sigma_,p,T)
    print("In projection_image_phi");
    print("beta_: ",beta_);
    print("sigma_:",sigma_);
    local phi_,nv,np,i:
        for np from 1 to T do 
            phi_[np][1]:=alpha[np]:
            for nv from 2 to num_var do 
                phi_[np][nv]:=beta_[nv-1]*alpha[np]-beta_[nv-1]*sigma_[1]+sigma_[nv] mod p:
            end do:
        end do:
    return [seq(convert(phi_[i],list),i=1..T)]:
end proc:

########################################################
# 3. GET_U - For vector interpolation
########################################################

get_u:=proc(M,col,alpha,p)
    print("In get_u");
    F:=[seq(convert(M[..,i],list),i=1..col)];
    U:=[seq(Interp(alpha,F[i],x) mod p,i=1..col)];
    return U;        
end proc:

########################################################
# 4. MQRFR - Modified QR for Rational Functions
########################################################

MQRFR:=proc(r0,r1,t0,t1,p)
    print("In MQRFR"):
    local r,t,q,i,f,g,qmax,lcg:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
    print("r0=",r0):
    # print("r1=",r1):
    # print("t0=",t0):
    # print("t1=",t1):
    f:=r0:
    g:=t1:
    qmax:=0:
    i:=1:
    while r[i] <> 0 do
        q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p:
        if degree(q[i],x)> qmax then 
            qmax:=degree(q[i],x):
            print("q[",i,"]=",q[i]," degree q=",qmax):
            f:=r[i]:
            g:=t[i]:
            # print("r[",i-1,"]=",r[i-1]):
            # print("q[",i,"]=",q[i]):
            # print("f=",f):
            # print("g=",g):
        end if:
        t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
        i:=i+1:
        if qmax <=1 or gcd(f,g) <> 1 or g = 0 then 
            FAIL:
        end if:
    end do:
    lcg:=lcoeff(g):
    # print("lcg=",lcg):
    return f/lcg mod p,g/lcg mod p,qmax,lcg :
end proc:

########################################################
# 5. NDSA - For Vector Black Boxes
########################################################

with(LinearAlgebra):

NDSA:=proc(B,sigma_,beta_,num_var,p,num_points)
    local correct_degree,T,alpha,m,phi_,_phi,Y,u,f,g,dq,lcg,np,nv,i,r,
          lin_sys,temp,result,count,M,row,col,DQ:
    print("In NDSA");
    correct_degree:=false:
    lin_sys:=false:
    T:=num_points:
    temp:=[]:
    result:=[]:
    count:=0:
    
    while(not(correct_degree)) do
        count:=count+1:
        print("T:= ",T):
        r:=rand(p):
        alpha:=[seq(r()+r() mod p,i=1..T)]:
        print("alpha: ",alpha):
        m:=expand(product(x-alpha[j],j=1..T)) mod p:
        print("m: ",m):
        _phi:=projection_image_phi(num_var,alpha,beta_,sigma_,p,T):
        Y:=[seq(B(_phi[i],p),i=1..T)]:
        M:=convert(Y,Matrix):
        # print("Matrix M: ",M):
        row,col:=Dimension(M):
        # print("row: ",row, " col: ",col):
        
        if row =1 then 
            u:=Interp(alpha,Y,x)mod p:
        else
            lin_sys:=true: 
            u:=get_u(M,col,alpha,p):
        end if:
        print("u: ",u):
        
        if lin_sys = false then  
            result:=[MQRFR(m,u,0,1,p)]:
            dq:=result[1][3]:
        else 
            for i from 1 to nops(u) do 
                temp:=[op(temp),MQRFR(m,u[i],0,1,p)]:
                result:=[op(result),temp]:
                temp:=[]:
            end do:
            print("result: ",result):
            DQ:=[seq(result[i][3],i=1..nops(result))]:
            print("DQ: ",DQ):
            dq:=min(DQ):
            print("Minimum dq: ",dq):
        end if:
        
        if dq > 1 then  
            print("Termination condition met"):    
            print("T:= ",T):
            print("num_points:= ",num_points):      
            if num_points <> T then  
                
                return result,T,lin_sys:
            else
                return result,lin_sys: 
            end if:
        else 
            print("MQRFR failed. Trying again with more points"):
            T:=T*2:
            result:=[]:
            DQ:=[]:
        end if:
        if(count= 4) then break: end if:
    end do:
end proc:

########################################################
# 6. BMEA - Berlekamp-Massey Extended Algorithm
########################################################

BMEA := proc(v::list,p::posint,Z::name)
    local n,m,R0,R1,V0,V1,i,Q:
    print("In BMEA");
    print("v=",v);    
    n := iquo( nops(v), 2 ):
    print("n=",n);
    m := 2*n-1:
    print("m=",m);
    R0 := Z^(2*n):
    print("R0=",R0);
    R1 := add( v[m+1-i]*Z^i, i=0..m ) mod p:
    print("R1=",R1);
    V0 := 0:
    V1 := 1:
    while n <= degree(R1,Z) do
        R0,R1 := R1,Rem(R0,R1,Z,'Q') mod p:
        print("R0=",R0);   
        print("R1=",R1);

        V0,V1 := V1,Expand(V0-Q*V1) mod p:
        print("V0=",V0);   
        print("V1=",V1);
    od:
    i := 1/lcoeff(V1,Z) mod p:
    return i*V1 mod p:
end:

########################################################
# 7. GENERATE MONOMIALS
########################################################

generate_monomials:=proc(roots_,num_var,prime_points,vars)
    local m,mm,i,j,counter,M_,rem: 
    M_:=Vector(numelems(roots_),0):
    print("In generate_monomials"):
    print("roots_=",roots_):
    
    for i from 1 to numelems(roots_) do
        if(roots_[i]=0)then 
            print("roots_[",i,"] = 0"):
            return FAIL: 
        end if:
        mm:=roots_[i]:
        m:=1:
        for j from 1 to numelems(prime_points) do
            counter:=0:
            while mm mod prime_points[j] = 0 do
                mm:=iquo(mm,prime_points[j],'rem'):
                counter:=counter+1:
            end do:
            m:=m*vars[j]^counter:
        end do:
        M_[i]:=m:
    end do:
    
    if mm<> 1 then 
        print("Warning: mm=",mm," (should be 1)"):
        return FAIL:
    end if:
    
    return convert(M_,list):
end proc:

########################################################
# 8. ZIPPEL VANDERMONDE SOLVER
########################################################

Zippel_Vandermonde_solver:=proc(y::list,terms::integer,roots_::list,lambda_,p::integer)
    local M,fin_coeff,q,q_lambda_inv,V_inv_b,i,j:
    print("In Zippel_Vandermonde_solver"):
    M:=lambda_ mod p:
    fin_coeff:=Vector(terms,0):
    for i from 1 to terms do
        q:=quo(M,Z-roots_[i],Z):
        q_lambda_inv:= 1/ Eval(q,Z=roots_[i]) mod p:
        V_inv_b:=0:
        for j from 1 to terms do
            V_inv_b:=V_inv_b+coeff(q,Z,j-1)*y[j] mod p:
        end do:
        fin_coeff[i]:=V_inv_b*q_lambda_inv mod p:
    end do:
    return convert(fin_coeff,list):
end proc:

########################################################
# 9. CONSTRUCT FINAL POLYNOMIAL
########################################################

construct_final_polynomial:=proc(coeff_,Monomials)
    print("In construct_final_polynomial"):
    local i,f:
    f:=0:
    for i from 1 to numelems(coeff_) do
        f:=f+coeff_[i]*Monomials[i]:
    end do:
    return f:
end proc:

########################################################
# 10. MAIN MRFI FOR VECTOR BLACK BOXES
########################################################

MRFI_Vector:=proc(B, num_vars::integer, num_eqn::integer, vars::list, p::integer)
    local i, j, k, Primes, shift_, sigma_, num_list, den_list, 
          numerator_done, denominator_done, T, T_old,
          mqrfr_results, lin_sys, num_points_mqrfr,
          F_current, G_current, deg_num, deg_den,
          lambda_num_list, lambda_den_list, terms_num_list, terms_den_list,
          R_num_list, R_den_list, Roots_num_list, Roots_den_list,
          num_mono_list, den_mono_list, coeff_num_list, coeff_den_list,
          final_num_list, final_den_list, r_, temp, common_den_flag,
          den_lc, den_lc_inv:
    
    print("========================================"):
    print("Starting MRFI_Vector"):
    print("Number of parameters:", num_vars):
    print("Number of equations:", num_eqn):
    print("========================================"):
    
    r_ := rand(p):
    Primes := [seq(ithprime(i), i=1..num_vars)]:
    shift_ := [seq(r_(), i=1..num_vars-1)]:
    sigma_ := []:
    
    # Initialize storage
    num_list := [seq([], i=1..num_eqn)]:
    den_list := [seq([], i=1..num_eqn)]:
    numerator_done := [seq(false, i=1..num_eqn)]:
    denominator_done := [seq(false, i=1..num_eqn)]:
    print("numerator_done: ",numerator_done):
    print("denominator_done: ",denominator_done):
    T := 4:
    T_old := 1:
    
    # Initial NDSA call
    mqrfr_results, num_points_mqrfr, lin_sys := NDSA(B, [seq(1, i=1..num_vars)], shift_, num_vars, p, T):
        print("numerator_done: ",numerator_done):
        print("denominator_done: ",denominator_done):
    if not lin_sys then
        error "Expected a linear system (vector output)":
    end if:
    
    # Process initial results
    F_current := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
    print("F_current: ",F_current):
    G_current := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
    print("G_current: ",G_current):

    for k from 1 to num_eqn do
        num_list[k] := [eval(F_current[k], x=1)]:
        
        den_list[k] := [eval(G_current[k], x=1)]:
    end do:
    print("num_list: ",num_list):
    print("den_list: ",den_list):
    deg_num := [seq(degree(F_current[i], x), i=1..nops(F_current))]:
    print("deg_num: ",deg_num):
    deg_den := [seq(degree(G_current[i], x), i=1..nops(G_current))]:
    print("deg_den: ",deg_den):
    num_points_mqrfr := max(op(deg_num)) + max(op(deg_den)) + 2:
    print("num_points_mqrfr: ",num_points_mqrfr):
    # Check for common denominator
    common_den_flag := true:
    print("G_current: ",G_current):
    for k from 2 to num_eqn do
        if G_current[k] <> G_current[1] then
            common_den_flag := false:
            break:
        end if:
    end do:
    print("Common denominator:", common_den_flag):
    
    # Main evaluation loop
    while true do
        print("T=", T):
        
        for j from T_old to 2*T-1 do
            sigma_ := [op(sigma_), [seq(Primes[i]^j mod p, i=1..nops(Primes))]]:
            
            # mqrfr_results, temp, lin_sys := NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr):
            mqrfr_results, lin_sys := NDSA(B, sigma_[j], shift_, num_vars, p, num_points_mqrfr):
            
            F_current := [seq(mqrfr_results[i][1], i=1..nops(mqrfr_results))]:
            print("F_current: ",F_current):
            G_current := [seq(mqrfr_results[i][2], i=1..nops(mqrfr_results))]:
            print("G_current: ",G_current):
            

            for k from 1 to num_eqn do
                print("k=",k):
                print("num_eqn=",num_eqn):
                print("numerator_done[",k,"]=",numerator_done[k]):
                print("num_list: ",num_list):
                print("den_list: ",den_list):
                if not numerator_done[k] then
                    num_list[k] := [op(num_list[k]), eval(F_current[k], x=sigma_[j][1]) mod p]:
                end if:
                if not denominator_done[k] then
                    den_list[k] := [op(den_list[k]), eval(G_current[k], x=sigma_[j][1]) mod p]:
                end if:
            end do:
        end do:
        print("num_list: ",num_list):
        
        # Apply BMEA
        lambda_num_list := []:
        terms_num_list := []:
        R_num_list := []:
        
        all_num_done := true:
        for k from 1 to num_eqn do
            print("k=",k):
            print("num_eqn=",num_eqn):
            print("numerator_done[",k,"]=",numerator_done[k]): 
            numerator_done[k] := false: 
            if not numerator_done[k] then
                temp := BMEA(num_list[k], p, Z):
                print("temp (numerator)=",temp):
                lambda_num_list := [op(lambda_num_list), temp]:
                print("lambda_num_list: ",lambda_num_list):
                terms_num_list := [op(terms_num_list), degree(temp, Z)]:
                print("terms_num_list: ",terms_num_list):
                R_num_list := [op(R_num_list), Roots(temp) mod p]:
                print("R_num_list: ",R_num_list):
                
                if nops(R_num_list[k]) > 0 and R_num_list[k][1][1] = 0 then
                    R_num_list[k] := remove(x->x=[0,1], R_num_list[k]):
                    terms_num_list[k] := terms_num_list[k] - 1:
                end if:
                
                if nops(R_num_list[k]) = terms_num_list[k] and terms_num_list[k] < T then
                    numerator_done[k] := true:
                else
                    all_num_done := false:
                end if:
            end if:
        end do:
        
        # Process denominators
        print("den_list: ",den_list):
        lambda_den_list := []:
        terms_den_list := []:
        R_den_list := []:
        
        all_den_done := true:
        print("common_den_flag: ",common_den_flag):
        if common_den_flag then
            print("denominator_done: ",denominator_done):
            denominator_done[1] := false:
            if not denominator_done[1] then
                temp := BMEA(den_list[1], p, Z):
                print("temp (denominator)=",temp):
                lambda_den_list := [temp]:
                print("lambda_den_list: ",lambda_den_list):
                terms_den_list := [degree(temp, Z)]:
                print("terms_den_list: ",terms_den_list):
                R_den_list := [Roots(temp) mod p]:
                print("R_den_list: ",R_den_list):
                
                if nops(R_den_list[1]) > 0 and R_den_list[1][1][1] = 0 then
                    R_den_list[1] := remove(x->x=[0,1], R_den_list[1]):
                    terms_den_list[1] := terms_den_list[1] - 1:
                end if:
                
                if nops(R_den_list[1]) = terms_den_list[1] and terms_den_list[1] < T then
                    for k from 1 to num_eqn do
                        denominator_done[k] := true:
                    end do:
                else
                    all_den_done := false:
                end if:
            end if:
        else
            for k from 1 to num_eqn do
                if not denominator_done[k] then
                    temp := BMEA(den_list[k], p, Z):
                    lambda_den_list := [op(lambda_den_list), temp]:
                    terms_den_list := [op(terms_den_list), degree(temp, Z)]:
                    R_den_list := [op(R_den_list), Roots(temp) mod p]:
                    
                    if nops(R_den_list[k]) > 0 and R_den_list[k][1][1] = 0 then
                        R_den_list[k] := remove(x->x=[0,1], R_den_list[k]):
                        terms_den_list[k] := terms_den_list[k] - 1:
                    end if:
                    
                    if nops(R_den_list[k]) = terms_den_list[k] and terms_den_list[k] < T then
                        denominator_done[k] := true:
                    else
                        all_den_done := false:
                    end if:
                end if:
            end do:
        end if:
        
        if all_num_done and all_den_done then
            print("All components recovered!"):
            break:
        end if:
        
        T_old := 2*T:
        T := T*2:
    end do:
    
    # Extract roots and build polynomials
    Roots_num_list := [seq([seq(r[1], r in R_num_list[k])], k=1..num_eqn)]:
    
    if common_den_flag then
        Roots_den_list := [[seq(r[1], r in R_den_list[1])]]:
        for k from 2 to num_eqn do
            Roots_den_list := [op(Roots_den_list), Roots_den_list[1]]:
        end do:
    else
        Roots_den_list := [seq([seq(r[1], r in R_den_list[k])], k=1..num_eqn)]:
    end if:
    
    # Generate monomials and coefficients
    num_mono_list := []:
    coeff_num_list := []:
    final_num_list := []:
    
    for k from 1 to num_eqn do
        temp := generate_monomials(Roots_num_list[k], num_vars, Primes, vars):
        if temp = FAIL then return FAIL: end if:
        num_mono_list := [op(num_mono_list), temp]:
        
        coeff_num_list := [op(coeff_num_list), 
            Zippel_Vandermonde_solver(num_list[k], terms_num_list[k], 
                                     Roots_num_list[k], lambda_num_list[k], p)]:
        
        final_num_list := [op(final_num_list), 
            construct_final_polynomial(coeff_num_list[k], num_mono_list[k])]:
    end do:
    
    # Generate denominators
    den_mono_list := []:
    coeff_den_list := []:
    final_den_list := []:
    
    if common_den_flag then
        temp := generate_monomials(Roots_den_list[1], num_vars, Primes, vars):
        if temp = FAIL then return FAIL: end if:
        
        coeff_den := Zippel_Vandermonde_solver(den_list[1], terms_den_list[1],
                                               Roots_den_list[1], lambda_den_list[1], p):
        
        temp_den := construct_final_polynomial(coeff_den, temp):
        final_den_list := [seq(temp_den, k=1..num_eqn)]:
    else
        for k from 1 to num_eqn do
            temp := generate_monomials(Roots_den_list[k], num_vars, Primes, vars):
            if temp = FAIL then return FAIL: end if:
            
            coeff_den := Zippel_Vandermonde_solver(den_list[k], terms_den_list[k],
                                                   Roots_den_list[k], lambda_den_list[k], p):
            
            final_den_list := [op(final_den_list), 
                construct_final_polynomial(coeff_den, temp)]:
        end do:
    end if:
    
    # Normalize
    for k from 1 to num_eqn do
        den_lc := lcoeff(final_den_list[k], order=grlex(seq(vars[i], i=1..num_vars))):
        den_lc_inv := 1/den_lc mod p:
        final_num_list[k] := final_num_list[k] * den_lc_inv mod p:
        final_den_list[k] := final_den_list[k] * den_lc_inv mod p:
    end do:
    
    print("========================================"):
    print("RECOVERY COMPLETE"):
    print("========================================"):
    
    return final_num_list, final_den_list:
end proc:

########################################################
# 11. DATA GENERATOR FOR TEST CASES
########################################################

get_data:=proc(test_case)
    local Sys,Vars,p,i,params;
    p:=2^31-1:
    if test_case = "bspline" then 
        Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1, 
                x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17, 
                -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
                -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2, 
                -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13, 
                -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
                2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21, 
                2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7, 
                2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9, 
                -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
                2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16, 
                2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:
    elif test_case ="small_sys_low_deg" then 
         Sys:={x1+y1*x2+y2*x3-1,y2*x1+x2+y1*x3-2,(y1-y2)*x1-x2+y2*x3-7}:
    elif test_case ="small_Sys" then
         Sys:={x1+y1*x2+y1-3,y2*x1+x2+y1-1}:
    elif test_case = "example" then
        # Your example: [x1 = -(y1^2-2*y1+3)/(y1*y2-1), x2 = -(y1*y2-y1-3*y2+1)/(y1*y2-1)]
        Sys := {(y1*y2-1)*x1 + (y1^2-2*y1+3), (y1*y2-1)*x2 + (y1*y2-y1-3*y2+1)}:
    end if:
    
    Vars := { seq( x||i, i=1..nops(Sys) )}:
    params := indets(Sys) minus Vars:
    return Sys,Vars,convert(params,list),p:
end proc:

print("Complete MRFI Vector Package Loaded Successfully!"):
print("========================================"):
