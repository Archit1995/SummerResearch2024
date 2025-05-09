Construct_Rational_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        local var,num,denom_,a,v;
        var:=vars:
        num:=f:
        denom_:=g:
        a:=num/denom_;
        return Eval(a,{seq(var[v]=point_[v],v=1..numelems(point_))}) mod p:
    end proc:
    return BB:
end proc: