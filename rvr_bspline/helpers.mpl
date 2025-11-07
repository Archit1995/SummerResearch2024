# Helper for solving the system
get_eqn:=proc(Sys,vars)option remember:
    print("in get_eqn"):
    return solve(Sys,Vars):
end proc:

########################################################
# 1. REORDERING HELPER FUNCTIONS
########################################################

reording:=proc(unordered_soln,num_eqn)
    local temp_var,component1,ordered_soln,i;
    ordered_soln:=[seq(0,i=1..num_eqn)];
    for i from 1 to num_eqn do 
        component1:=get_component(op(1,unordered_soln[i])):
        ordered_soln[component1]:=unordered_soln[i];
    end do;
    return ordered_soln;
end proc:

get_component:=proc(expression)
    local temp_var;
    temp_var:=convert(expression,string);
    return parse(temp_var[2..length(temp_var)]);
end proc:

########################################################
# 3. GET_U - For vector interpolation
########################################################

get_u:=proc(M,col,alpha,p)
    # print("In get_u");
    local F,U,i:
    F:=[seq(convert(M[..,i],list),i=1..col)];
    # print("F: ",F);
    U:=[seq(Interp(alpha,F[i],x) mod p,i=1..col)];
    # F:=convert(M[..,col][..numelems(alpha)],list);
    # U:=[seq(Interp(alpha,x) mod p)];
    return U;        
end proc:

Deterministic_get_u:=proc(M,col,alpha,p)
    # print("In get_u deterministic"):
    local F,i,U:
    F:=convert(M[..,col][..numelems(alpha)],list);
    # print("F: ",F);
    U:=Interp(alpha,F,x) mod p;
    return U;        
end proc:
