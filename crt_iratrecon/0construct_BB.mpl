Construct_Rational_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        # print("point_=",point_):
        local var,num,denom_,a,v;
        var:=vars:
        num:=f:
        denom_:=g:
        a:=num/denom_;
        if Eval(denom_,{seq(var[v]=point_[v],v=1..numelems(point_))}) mod p = 0 then
            # denominator is zero
            error "Denominator is zero":
        else 
            return Eval(a,{seq(var[v]=point_[v],v=1..numelems(point_))}) mod p:
        end if:
    end proc:
    return BB:
end proc: