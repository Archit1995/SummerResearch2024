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
p:=2^31-1;
vars:={x,y,z};
f:=randpoly(vars,sparse,degree=20):
B:=Construct_Blackbox(f,vars);
print(B);
whattype(B);