# 1. Have a black box for f/g
Construct_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        local u,v,var,num,den,a;
        var:=vars:
        # num:=f:
        # denom_:=g:
        # a:=num/denom_:
        a:=f:
        # print("a: ",a);
        return modp([seq(eval(a,{seq(var[v]=point_[u][v],v=1..numelems(point_[u]))}),u=1..numelems(point_))],p):
    end proc:
    return BB:
end proc:
get_degree_by_variable:=proc(vars,B)
    local num_var,temp,eval_point,deg,i,j;
    num_var:=numelems(vars);
    print("vars = ",vars);
    if num_var =1 then 
        temp:=B([[vars[1]]],p);
        print("temp = ",temp[1]);
        return degree(temp[1],vars[1]);
    end if;   
    for i from 1 to num_var do
        eval_point:=convert(Vector(num_var),list);
        for j from 1 to num_var do 
            # pos:=j mod num_var
            if i =j then 
                eval_point[j]:=vars[j]; 
            else 
                eval_point[j]:=j;
            end if;
        end do;
        print("eval point = ",eval_point);
        temp:=B([eval_point],p);
        deg[i]:=degree(temp[1],vars[i]);
        print(temp[1]);
    end do;
    return convert(deg,list);
end proc:
p:=19;
num_var:=2:
# num_var:=1;
vars:={x,y}:
# vars:={x};
# prime_points:=generate_evaulation_primes(num_var):
# prime_points:
# p:=2^31-1:
# p:=19;
# f:=randpoly(vars,sparse,degree=7) mod p:
# g:=randpoly(vars,sparse,degree=5) mod p:
# B:=Construct_Blackbox(f,g,vars);
# print(B);
# get_degree_by_variable(vars,B);


vars:={w,x,y,z};

f:=randpoly(vars,sparse,degree=70) mod p;
g:=randpoly(vars,sparse,degree=5) mod p:
B:=Construct_Blackbox(f,g,vars);
get_degree_by_variable(vars,B);