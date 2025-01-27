Construct_Blackbox:=proc(f,vars)
    local BB:
    BB:=proc(point_,p)
        local var,a;
        var:=vars:
        a:=f:
        # print("a: ",a);
        return modp([seq(eval(a,{seq(var[v]=point_[u][v],v=1..numelems(point_[u]))}),u=1..numelems(point_))],p):
    end proc:
    return BB:
end proc:

# get_degree_by_variable:=proc(vars,B)
#     local num_var,temp,eval_point,deg,i,j;
#     num_var:=numelems(vars);
#     print("vars = ",vars);
#     if num_var =1 then 
#         temp:=B([[vars[1]]],p);
#         print("temp = ",temp[1]);
#         return degree(temp[1],vars[1]);
#     end if;
#     for i from 1 to num_var do
#         eval_point:=convert(Vector(num_var),list);
#         for j from 1 to num_var do 
#             if i =j then 
#                 eval_point[j]:=vars[j]; 
#             else 
#                 eval_point[j]:=j;
#             end if;
#         end do;
#         print("eval point = ",eval_point);
#         temp:=B([eval_point],p);
#         deg[i]:=degree(temp[1],vars[i]);
#         # print(temp[1]);
#     end do;
#     return convert(deg,list);
# end proc:

eea_gcd:=proc(r0,r1,t0,t1,p)
    local r,t,q,i,f,g,qmax:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
     print("t0= ",t0):
     print("t1= ",t1):
    f:=r0:
    g:=t1:
    qmax:=1:
    i:=1:
    while r[i] <> 0 do
    #  1. find quotient and remainder
    q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p:
     print("r[",i,"]=",r[i]):
     print("q[",i,"]=",q[i]):
     print("t[",i,"]=",t[i]):
     print("Degree of q[",i,"]=",degree(q[i],x)):
     print("--------------------------------------"):
     print("r[",i+1,"]=",r[i+1]):
    if degree(q[i],x)>= qmax then 
        qmax:=degree(q[i],x):
        f:=r[i]:
        g:=t[i]:
    end if:
    #   reduction step to update t.
    t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
    #  Normalization step for r and t 
     if r[i+1]<> 0 then 
        rho:=lcoeff(r[i+1]):
        r[i+1]:=r[i+1]/rho mod p:
        t[i+1]:=t[i+1]/rho mod p:
    end if:
     print("f=",f):
     print("g=",g):
     print("______________________________________"):
    i:=i+1:
    end do:
    return f,g:
end proc:
