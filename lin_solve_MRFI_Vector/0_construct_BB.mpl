########################################################
# 0. CONSTRUCT BLACK BOX - For Linear Systems
########################################################

Constuct_Sys_Blackbox:=proc(Sys,Vars,params) 
    local Lin_BB:
    Lin_BB:=proc(point_,p)option remember:
        local unordered_soln,Soln,soln,var,i,v:
        # print("point_",point_):
        global counter;
        counter:=counter+1;
        var:=params:
        Soln:=get_eqn(Sys,Vars):
        unordered_soln:=convert(Soln,list):
        soln:=reording(unordered_soln,nops(Sys)):
        # print("soln=",soln);
        return [seq(eval(op(2,soln[i]),{seq(var[v]=point_[v],v=1..numelems(point_))}),i=1..nops(Sys))] mod p:
    end proc:
    return Lin_BB:
end proc:

Construct_Rational_Blackbox:=proc(f,g,vars)
    local BB:
    BB:=proc(point_,p)
        local var,num,denom_,a,v;
        global counter;
        counter:=counter+1;       
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



