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
end proc;