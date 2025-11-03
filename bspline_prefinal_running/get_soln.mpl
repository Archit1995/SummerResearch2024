
read "./data_gen.mpl":
read "./helpers.mpl":
Sys, Vars, params, p := get_data("bspline"):
get_soln:=proc(Sys,Vars,params)
    local unordered_soln,Soln,soln,var,i,v:
    var:=params:
    Soln:=get_eqn(Sys,Vars):
    unordered_soln:=convert(Soln,list):
    soln:=reording(unordered_soln,nops(Sys)):
    return soln:
end proc:

x:=get_soln(Sys,Vars,params);
for i from 1 to nops(x) do
    lprint(" x[",i,"]=: ",x[i]);
end do;