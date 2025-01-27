# 1. Have a black box for f/g. 
Construct_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        # print("point_ ",point_);
        local var,num,denom_,a,i:
        var:=vars:
        num:=f:
        denom_:=g:
        a:=num/denom_:
        return Eval(a,{seq(var[i]=point_[i],i=1..numelems(point_))}) mod p: 
    end proc:
    return BB:
end proc:
generate_random_vector:=proc(n,p)
    local r,i:
    r:=rand(p):
    return [seq(r(),i=1..n)]:
end proc:

# MQRFR(m,u,0,1,p)
MQRFR:=proc(r0,r1,t0,t1,p)
    local r,t,q,i,f,g,qmax,lcg:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
    f:=r0:
    g:=t1:
    qmax:=0:
    i:=1:
    while r[i] <> 0 do
    #  1. find quotient and remainder
    q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p:
    if degree(q[i],x)>= qmax then 
        qmax:=degree(q[i],x):
        print("max degree = ",qmax);
        f:=r[i]:
        g:=t[i]:
    end if:
    #   reduction step to update t.
    t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
    i:=i+1:
    if qmax <=1 or gcd(f,g) <> 1 then 
        FAIL:
    end if:
    end do:
    lcg:=lcoeff(g):
    return f/lcg mod p,g/lcg mod p,qmax,lcg:# make g monic. Also add a check to see that g !=0
end proc:

Early_termination_seperation:=proc(B,p)
    local r,u,f,g,dq,num_points,correct_degree,points,Y,m,check,lc_u,sigma_,shift_,alpha,phi_,nv,np,m_u,i,lcg:
    r:=rand(p):
    num_points:=1:
    correct_degree:=true:
    while(correct_degree) do
        # print("num_points: ",num_points):
        r:=rand(p):
        alpha:=[seq(r(),i=1..num_points)]:
        sigma_:=generate_random_vector(numelems(vars),p):
        # print("sigma_ ",sigma_);
        shift_:=generate_random_vector(numelems(vars)-1,p):
        # print("shift_ = ",shift_);
        temp1:=Vector(num_var);
        temp2:=convert(temp1,list);
        phi_:=[seq(temp2,i=1..num_points)]:
        # print("phi_",phi_);
        for np from 1 to num_points do 
            phi_[np][1]:=alpha[np]:
            for nv from 2 to num_var do 
                phi_[np][nv]:=shift_[nv-1]*alpha[np]-shift_[nv-1]*sigma_[1]+sigma_[2] mod p:
            end do:
        end do:
        # print("phi_ ",phi_);

        Y:=[seq(B(phi_[i],p),i=1..num_points)]:
        m:=Expand(product(x-alpha[i],i=1..num_points)) mod p:
        print("m",m);
        u:=Interp(alpha,Y,x)mod p:
        print("u =",u);
        m_u:=u/lc_u mod p:
        f,g,dq,lcg:=MQRFR(m,u,0,1,p):
        if dq > 1 then 
            break:
        else 
            num_points:=num_points*2:
        end if:
    end do:
    return f,g,lcg,sigma_,shift_,num_points:
end proc:

Monte_Carlo_seperation:=proc(num_deg,denom_deg,vars,p)
    local ff,gg,B,f,g,lc_g,sigma_,shift_,Phi,f_x,g_x,np:
    ff:=randpoly(vars,sparse,degree=num_deg) mod p:
    gg:=randpoly(vars,sparse,degree=denom_deg) mod p:
    gg:=gg/lcoeff(gg) mod p:
    print("ff ",ff):
    print("gg ",gg):
    B:=Construct_Blackbox(ff,gg,vars):
    f,g,lc_g,sigma_,shift_,np:=Early_termination_seperation(B,p):
    print("f = ",f);
    print("g = ",g);
    # Phi:=shift_[1]*x-shift_[1]*sigma_[1]+sigma_[2] mod p:j=2 ..num_var
    Phi:=[seq(modp(shift_[j-1]*x-shift_[j-1]*sigma_[1]+sigma_[j-1],p),j=2..num_var)];
    print("Phi ",Phi):
    f_x:=Expand(subs(vars[2]=Phi[1],ff)) mod p:
    g_x:=Expand(subs(vars[2]=Phi[1],gg)) mod p:
    f_x:=f_x/lcoeff(g_x) mod p:
    g_x:=g_x/lcoeff(g_x) mod p: 
    print("f_x =",f_x);
    print("g_x =",g_x);
    return f_x-f,g_x-g,f,g,np:
    # print("____________________________________________________________");
end proc:
vars:={x,y}:
# vars:={x1,x2}:
num_var:=numelems(vars):
p:=2^31-1:
# Test 1
tester:=proc(num_deg,denom_deg,vars,p)
    local diff_num,diff_denom,counter,num_recon,denom_recon,num_Points:
    counter:=0:
    while true do 
        diff_num,diff_denom,num_recon,denom_recon,num_Points:=Monte_Carlo_seperation(num_deg,denom_deg,vars,p);
        counter:=counter+1:
        break;
        # print("diff_num: ",diff_num):
        # print("diff_denom: ",diff_denom):
        if diff_num = 0 and diff_denom = 0 then 
            return counter,num_recon,denom_recon,num_Points:

        end if:
    end do:
end proc:
# Test 1
# num_deg:=20:
# denom_deg:=10:
# tester(num_deg,denom_deg,vars,p);
# Test 2
# num_deg:=35:
# denom_deg:=15:
# tester(num_deg,denom_deg,vars,p);
# # Test 3
# num_deg:=44:
# denom_deg:=30:
# tester(num_deg,denom_deg,vars,p);
# # Test 4
# num_deg:=80:
# denom_deg:=40:
# tester(num_deg,denom_deg,vars,p);
# Test 4
# num_deg:=244:
num_deg:=9:
# denom_deg:=130:
denom_deg:=4:
count_,NUM,DENOM,NP:=tester(num_deg,denom_deg,vars,p);
print("count_: ",count_);
print("num points used: ",NP);
print("NUM: ",NUM);
print("DENOM: ",DENOM);
