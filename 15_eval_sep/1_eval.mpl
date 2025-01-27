# The key idea for this is to check that 
Construct_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        # #print("point_ ",point_);
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
MQRFR:=proc(r0,r1,t0,t1,p,vars)
    local r,t_,q,i,f,g,qmax,lcg:
    r[0]:=r0:
    r[1]:=r1:
    t_[0]:=t0:
    t_[1]:=t1:
    f:=r0:
    g:=t1:
    qmax:=0:
    i:=1:
    while r[i] <> 0 do
    #  1. find quotient and remainder
    q[i]:= Quo(r[i-1],r[i],vars[1],'r[i+1]') mod p:
    print("r[i+1] = ",r[i+1]);
    if degree(q[i],vars[1])>= qmax then 
        qmax:=degree(q[i],vars[1]):
        # print("max degree = ",qmax);
        f:=r[i]:
        g:=t_[i]:
        print("g = ",g);
    end if:
    #   reduction step to update t.
    t_[i+1]:=Expand(t_[i-1]-q[i]*t_[i])mod p:

    i:=i+1:
    if qmax <=1 or gcd(f,g) <> 1 then 
        FAIL:
    end if:
    end do:
    lcg:=lcoeff(g,vars[1]):
    # print("lcg = ",lcg);
    return f/lcg mod p,g/lcg mod p,qmax,lcg:# make g monic. Also add a check to see that g !=0
end proc:

Early_termination_seperation:=proc(B,p,vars)
    local r,u_,f,g,dq,num_points,correct_degree,points,Y,m,check,lc_u,sigma_,shift_,alpha,phi_,nv,np,m_u,i,lcg:
    r:=rand(p):
    num_points:=1:
    correct_degree:=true:
    while(correct_degree) do
        print("num_points: ",num_points):
        r:=rand(p):
        alpha:=[seq(r(),i=1..num_points)]:
        print("alpha = ",alpha);
        sigma_:=generate_random_vector(numelems(vars),p):
        print("sigma_ ",sigma_);
        shift_:=generate_random_vector(numelems(vars)-1,p):
        print("shift_ = ",shift_);
        temp1:=Vector(num_var);
        temp2:=convert(temp1,list);
        temp3:=[seq(temp2,i=1..num_points)]:
        phi_:=convert(temp3,Array):
        # #print("phi_",phi_);
        for np from 1 to num_points do 
            # phi_[np][1]:=alpha[np]:
            phi_[np,1]:=alpha[np]:
            for nv from 2 to num_var do 
            # for nv from 1 to num_var-1 do 
            # Phi:=[seq(modp(shift_[j]*x-shift_[j]*sigma_[1]+sigma_[j+1],p),j=1..num_var-1)];
                phi_[np,nv]:=shift_[nv-1]*alpha[np]-shift_[nv-1]*sigma_[1]+sigma_[nv] mod p:
                # phi_[np,nv]:=shift_[nv]*alpha[np]-shift_[nv]*sigma_[1]+sigma_[nv+1] mod p:
            end do:
        end do:
        print("phi_ ",phi_);

        Y:=[seq(B(phi_[i],p),i=1..num_points)]:
        m:=Expand(product(vars[1]-alpha[i],i=1..num_points)) mod p:
        # m:=Expand(product(x-alpha[i],i=1..num_points)) mod p:
        print("m",m);
        u_:=Interp(alpha,Y,vars[1])mod p:
        # u_:=Interp(alpha,Y,x)mod p:
        print("u =",u_);
        # m_u:=u_/lcoeff(u_) mod p:
        f,g,dq,lcg:=MQRFR(m,u_,0,1,p,vars):
        # f,g,dq,lcg:=MQRFR(m,m_u,0,1,p,vars):
        # if(num_points >32) then break; end if;
        if dq > 1 then 
            break:
        else 
            num_points:=num_points*2:
        end if:
        print("____________________________________________________________");
    end do:
    return f,g,lcg,sigma_,shift_,num_points,phi_,alpha:
end proc:

Monte_Carlo_seperation:=proc(num_deg,denom_deg,vars,p)
    local i,j,ff,gg,B,f,g,lc_g,sigma_,shift_,Phi,f_x,g_x,np:
    ff:=randpoly(vars,sparse,degree=num_deg) mod p:
    gg:=randpoly(vars,sparse,degree=denom_deg) mod p:
    gg:=gg/lcoeff(gg) mod p:
    # ff:=ff/lcoeff(gg) mod p:
    #print("ff ",ff):
    #print("gg ",gg):
    B:=Construct_Blackbox(ff,gg,vars):
    f,g,lc_g,sigma_,shift_,np,phi_,alpha:=Early_termination_seperation(B,p,vars):
    Phi:=[seq(modp(shift_[j]*vars[1]-shift_[j]*sigma_[1]+sigma_[j+1],p),j=1..num_var-1)];# IS just for checking the correctness of the code.
    #print("f = ",f);
    #print("g = ",g);
    # Phi:=shift_[1]*x-shift_[1]*sigma_[1]+sigma_[2] mod p:j=2 ..num_var- 
    print("===============================================================");
    print("Phi ",Phi):
    test_Phi_alpha:=[seq(modp(eval(Phi[1],vars[1]=alpha[i]),p),i=1..numelems(alpha))]:
    print("checking Phi(alpha)",test_Phi_alpha);
    f_x:=Expand(eval(ff,[seq(vars[i+1]=Phi[i],i=1..num_var-1)])) mod p:
    g_x:=Expand(eval(gg,[seq(vars[i+1]=Phi[i],i=1..num_var-1)])) mod p:
    certificate:=lcoeff(g_x,vars[1]);
    print("lcoeff(g_x)= ",certificate/shift_[1]mod p);
    print("lg_c = ",lc_g/shift_[1]mod p);
    # print("lcoeff(g_x)/lg_c= ",lcoeff(g_x,vars[1])/lc_g mod p);    
    f_x:=f_x/lcoeff(g_x) mod p:
    g_x:=g_x/lcoeff(g_x) mod p: 
    print("f_x =",f_x);
    print("g_x =",g_x);
    print("===============================================================");
    return f_x-f,g_x-g,f,g,np,ff,gg,lc_g,phi_,sigma_,certificate:
    
end proc:
vars:={x,y}:
# vars:={x,y}:
# vars:={x1,x2}:
num_var:=numelems(vars):
# p:=2^31-1:
p:=509:
# Test 1
tester:=proc(num_deg,denom_deg,vars,p)
    local diff_num,diff_denom,counter,num_recon,denom_recon,num_Points:
    counter:=0:
    while true do 
        diff_num,diff_denom,num_recon,denom_recon,num_Points,ff,gg,const_,phi_,sigma_,certificate:=Monte_Carlo_seperation(num_deg,denom_deg,vars,p);
        counter:=counter+1:
        # break;
        # #print("diff_num: ",diff_num):
        #print("diff_denom: ",diff_denom):
        # eval_A:=Eval(ff,[seq(vars[i]=sigma_[i],i=1..num_var)]) mod p:
        # eval_B:=Eval(gg,[seq(vars[i]=sigma_[i],i=1..num_var)]) mod p:
        # print("f(sigma1,sigma2,sigma3) = ",eval_A);
        # print("g(sigma1,sigma2,sigma3) = ",eval_B);
        # ev_n1:=Eval(NUM,vars[1]=sigma_[1]^exp_[1]) mod p;
        # ev_d1:=Eval(DENOM,vars[1]=sigma_[1]^exp_[1]) mod p;
        # ev_n1:=Eval(num_recon,vars[1]=sigma_[1]) mod p;
        # ev_d1:=Eval(denom_recon,vars[1]=sigma_[1]) mod p;
        if diff_num = 0 and diff_denom = 0 then  
            return counter,num_recon,denom_recon,num_Points,ff,gg,const_,phi_,sigma_:
            

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
num_deg:=2:
# denom_deg:=130:
denom_deg:=3:
count_,NUM,DENOM,NP,ff,gg,lc_denom,phi_,sigma_:=tester(num_deg,denom_deg,vars,p);
print("count_: ",count_);
# print("num points used: ",NP);
# print("NUM: ",NUM);
# print("DENOM: ",DENOM);
# print("A: ",ff);
# print("B_: ",gg);
# print("lc_denom: ",lc_denom);
# lc_denom_inv:=1/lc_denom mod p;
# print("lc_denom_inv: ",lc_denom_inv);
# print("sigma_: ",sigma_);
# print("phi_: ",phi_);
# # exp_:=[1,2,3]:
# # eval_ff:=Eval(ff,[seq(vars[i]=sigma_[i]^exp_[i],i=1..num_var)]) mod p:
# # eval_gg:=Eval(gg,[seq(vars[i]=sigma_[i]^exp_[i],i=1..num_var)]) mod p:


for verifier from 1 to 20 do 
    anchor:=generate_random_vector(numelems(vars),p):
    eval_ff:=Eval(ff,[seq(vars[i]=anchor[i],i=1..num_var)]) mod p:
    eval_gg:=Eval(gg,[seq(vars[i]=anchor[i],i=1..num_var)]) mod p:
    ev_n1:=Eval(NUM,vars[1]=anchor[1]) mod p;
    ev_d1:=Eval(DENOM,vars[1]=anchor[1]) mod p;
    print("verifier: ",verifier);
    print(eval_ff/ev_n1 mod p);
    print(eval_gg/ev_d1 mod p);
    print("____________________________________________________________");
end do:

# print("f(sigma1,sigma2,sigma3) = ",eval_A);
# print("g(sigma1,sigma2,sigma3) = ",eval_B);
# # ev_n1:=Eval(NUM,vars[1]=sigma_[1]^exp_[1]) mod p;
# # ev_d1:=Eval(DENOM,vars[1]=sigma_[1]^exp_[1]) mod p;
# ev_n1:=Eval(NUM,vars[1]=sigma_[1]) mod p;
# ev_d1:=Eval(DENOM,vars[1]=sigma_[1]) mod p;
# print("NUM(sigma1) = ",ev_n1);
# print("DENOM(sigma1) = ",ev_d1);
# print(eval_ff/ev_n1 mod p);
# print(eval_gg/ev_d1 mod p);