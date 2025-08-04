

Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1, x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17, -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
 -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2, -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13, -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
  2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21, 2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7, 2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9, -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
   2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16, 2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:

# Sys:={x1+y1*x2+y2^2*x3-1,y2^3*x1+x2+y1*x3-2,x1-(y1^2-y2)*x2+y2^3*x3-7};

Vars := { seq(   x||i, i=1..nops(Sys) )}:
params := indets(Sys) minus Vars:


Soln:=solve(Sys,Vars):
unordered_soln:=convert(Soln,list);

seq(whattype(unordered_soln[i]),i=1..nops(Sys));
type(unordered_soln[1],equation);


# temp_pos:=convert(component,integer);
# op(2,expand(soln[3]));
# op(2,(soln[3]));
get_component:=proc(expression)
    local temp_var;
    temp_var:=convert(expression,string);
    return parse(temp_var[2..length(temp_var)]);
end proc;

reording:=proc(unordered_soln,num_eqn)
    local temp_var,component1,ordered_soln;
    ordered_soln:=[seq(0,i=1..num_eqn)];
    for i from 1 to num_eqn do 
        component1:=get_component(op(1,unordered_soln[i])):
        ordered_soln[component1]:=unordered_soln[i];
    end do;
    return ordered_soln;
end proc:



soln:=reording(unordered_soln,nops(Sys));




# 1/f(x)=f(x)^-1
# f(x)/g(x)=f(x)*(g(x)^-1)
get_num_den:=proc(soln,num_eqn)
    local num,den,n;
    num_soln:=[];
    den_soln:=[];
    n=num_eqn;
    for i from 1 to num_eqn do 
        expression:=op(2,(soln[i]));
        print("nops(expression) = ",nops(expression));
        print("soln[",i,"]= ",expression):
        if type(expression,`^`) then 
            num_soln:=[op(num_soln),1];
            den_soln:=[op(den_soln),op(1,expression)];
        end if;
        if type(expression,`*`)then 
            num_soln:=[op(num_soln),op(1,expression)];
            den_soln:=[op(den_soln),op(1,op(2,expression))];
        end if;
    end do;
    return num_soln,den_soln;
end proc:

NUM_soln,DEN_soln:=get_num_den(soln,nops(Sys));
NUM_soln[14];
DEN_soln[14];
whattype(op(2,soln[14]));
nops(op(2,soln[14]));
op(1,op(2,soln[14]));
op(2,op(2,soln[14]));
op(3,op(2,soln[14]));
expand(soln[14]);
p:=19;
Eval(soln[14],{y1=1,y2=2,y3=3,y4=4,y5=5,y6=6}) mod p;
ep:=[seq(eval(op(2,soln[i]),{y1=1,y2=2,y3=3,y4=4,y5=5,y6=6}),i=1..nops(Sys))];
op(2,ep[1]);