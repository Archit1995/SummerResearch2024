get_data:=proc(test_case)
    local Sys,Vars,p,i,params;
    p:=2^31-1:
    if test_case = "bspline" then 
        Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1, x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17, -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
            -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2, -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13, -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
            2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21, 2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7, 2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9, -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
            2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16, 2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:
    elif test_case ="small_sys_low_deg" then 
         Sys:={x1+y1*x2+y2*x3-1,y2*x1+x2+y1*x3-2,(y1-y2)*x1-x2+y2*x3-7}:

    
    elif test_case ="small_Sys" then
         Sys:={x1+y1*x2+y2^2*x3-1,y2^3*x1+x2+y1*x3-2,x1-x2+y2^3*x3-7}:
        
    end if:
    Vars := { seq(   x||i, i=1..nops(Sys) )}:
    params := indets(Sys) minus Vars:
    return Sys,Vars,convert(params,list),p:
end proc: