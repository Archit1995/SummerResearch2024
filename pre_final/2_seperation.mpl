generate_random_vector:=proc(n,p)# generate anchor points p_1..p_n and shift B_2..B_n, may be used 
# to create  random exponents of anchor points. 
    local r,i:
    r:=rand(p):
    return [seq(r(),i=1..n)]:
end proc:
# Monte Carlo algorithm for early termination separation
Early_termination_seperation:=proc(B,p)
    local r,u,f,g,dq,num_points,correct_degree,points,Y,m,check,lc_u,anchor_points,shift_,alpha,phi_,nv,np,m_u,i,lcg:
    # r:=rand(p):
    num_points:=1:
    correct_degree:=true:
    # anchor_points:=generate_random_vector(numelems(vars),p):
    # print("anchor_points: ",anchor_points):
    shift_:=generate_random_vector(numelems(vars)-1,p):
    print("shift_: ",shift_):
    while(correct_degree) do
        print("num_points: ",num_points):
        r:=rand(p):
        alpha:=[seq(r(),i=1..num_points)]:
        print("alpha: ",alpha):
        anchor_points:=generate_random_vector(numelems(vars),p):
        print("anchor_points: ",anchor_points):
        print("f(var[1],shift[1]*var[1]-shift[1]*anchor_points[1]+anchor_points[2]"): 
        # _phi:=[seq([seq(0,nv=1..num_var)],np=1..num_points)]:\
        # phi_:=convert(_phi,Array):
        # print("phi_: ",phi_):
        # phi_[np][1]=alpha[np]:
        # f(x,y)=T(x,B2x-B2p1+p2)= T(x,Bixi-Bip1+p2)
        # f(p1,p2)=T(p1)
        for np from 1 to num_points do # projection ring_homomorphism
            phi_[np][1]:=alpha[np]:
            for nv from 2 to num_var do 
                phi_[np][2]:=shift_[nv-1]*alpha[np]-shift_[nv-1]*anchor_points[1]^np+anchor_points[2]^np mod p:
            end do:
        end do:
        _phi:=[seq(convert(phi_[i],list),i=1..num_points)]:
        print("_phi: ",_phi):
        # print("phi_[1]: ",):
        
        Y:=[seq(B(_phi[i],p),i=1..num_points)]:
        # if(num_points=16)then break: end if:
        # Y:=[B(phi_[1],p)]:
        print("Y: ",Y):
        Phi:=shift_[1]*x-shift_[1]*anchor_points[1]+anchor_points[2] mod p:
        print("Phi: ",Phi):
        m:=Expand(product(x-alpha[i],i=1..num_points)) mod p:
        print("m: ",m):
        u:=Interp(alpha,Y,x)mod p:
        print("u: ",u):
        # checking:=[seq(Eval(u,x=alpha[i]) mod p,i=1..num_points)]:
        # print("Y: ",Y):
        # print("checking: ",checking):
        m_u:=u/lc_u mod p:
        f,g,dq,lcg:=MQRFR(m,u,0,1,p):
        # print("In early termination seperation"):
        # print("f= ",f):
        # print("g= ",g):
        # print("lcg= ",lcg):
        if num_points=16 then break: end if:
        if dq > 1 then 
            break:
        else 
            num_points:=num_points*2:
        end if:
        print("====================================================="):
    end do:
    return f,g,lcg,anchor_points,shift_:
end proc:

test_separation:=proc(ff,gg,B,p)
    # ff:=randpoly(vars,sparse,degree=num_deg) mod p:
    # gg:=randpoly(vars,sparse,degree=denom_deg) mod p:
    # B:=Construct_Blackbox(ff,gg,vars):
    f,g,lc_g,anchor_points,shift_:=Early_termination_seperation(B,p):
    Phi:=shift_[1]*x-shift_[1]*anchor_points[1]+anchor_points[2] mod p:
    f_x:=Expand(subs(y=Phi,ff)) mod p:
    g_x:=Expand(subs(y=Phi,gg)) mod p:
    lc_g_x:=lcoeff(g_x):
    f_x:=f_x/lc_g_x mod p:
    g_x:=g_x/lc_g_x mod p:
    print("leading coefficient of g: ",lc_g):
    print("1/lc_g = ",1/lc_g mod p):
    print("leading coefficient of g_x: ",lc_g_x):
    print("lc_g*lc_g_x =    ",lc_g*lc_g_x mod p):
    print("1/lc_g_x = ",1/lc_g_x mod p):
    print("checking numerator: ",f-f_x):
    print("checking denominator: ",g-g_x):
    B_verify:=Construct_Blackbox(f,g,vars):
    print("B_verify: ",B_verify):
    
end proc: