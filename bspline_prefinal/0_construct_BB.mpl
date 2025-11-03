########################################################
# 0. CONSTRUCT BLACK BOX - For Linear Systems
########################################################

Constuct_Sys_Blackbox:=proc(Sys,Vars,params) 
    local Lin_BB:
    Lin_BB:=proc(point_,p)option remember:
        local unordered_soln,Soln,soln,var,i,v:
        # print("point_",point_):
        var:=params:
        Soln:=get_eqn(Sys,Vars):
        unordered_soln:=convert(Soln,list):
        soln:=reording(unordered_soln,nops(Sys)):
        # print("soln=",soln);
        return [seq(eval(op(2,soln[i]),{seq(var[v]=point_[v],v=1..numelems(point_))}),i=1..nops(Sys))] mod p:
    end proc:
    return Lin_BB:
end proc:

Constuct_Vector_Blackbox:=proc(vector_,Vars,params) 
    local VBB:
    VBB:=proc(point_,p)option remember:
    var:=params:
        return [seq(eval(op(2,vector_[i]),{seq(var[v]=point_[v],v=1..numelems(point_))}),i=1..nops(vector_))] mod p:
    end proc:
    return VBB:
end proc:

