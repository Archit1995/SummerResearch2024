
########################################################
# 11. DATA GENERATOR FOR TEST CASES
########################################################

get_data:=proc(test_case)
    print("in get_data"):
    local Sys,Vars,i,params,ff,gg;
    if nargs = 1 then
        if test_case = "bspline" then 
            Sys := {x7 + x12 - 1, x8 + x13 - 1, x21 + x6 + x11 - 1, 
                    x1*y1 + x1 - x2, x11*y3 + x11 - x12, x16*y5 - x17*y5 - x17, 
                    -x20*y3 + x21*y3 + x21, x3*y2 + x3 - x4,
                    -x8*y4 + x9*y3 + x9, 2*x1*y1^2 - 2*x1 - 2*x10 + 4*x2, 
                    -x10*y2 + x18*y2 + x18 - x19, 2*x11*y3^2 - 2*x11 + 4*x12 - 2*x13, 
                    -x13*y4 + x14*y4 + x14 - x15, 2*x15*y5^2 - 4*x16*y5^2 + 2*x17*y5^2 - 2*x17,
                    2*x19*y3^2 - 4*x20*y3^2 + 2*x21*y3^2 - 2*x21, 
                    2*x3*y2^2 - 2*x3 + 4*x4 - 2*x5, -x5*y3 + x6*y3 + x6 - x7, 
                    2*x7*y4^2 - 4*x8*y4^2 + 2*x9*y4^2 - 2*x9, 
                    -4*x10*y2^2 + 2*x18*y2^2 + 2*x2*y2^2 - 2*x18 + 4*x19 - 2*x20,
                    2*x12*y4^2 - 4*x13*y4^2 + 2*x14*y4^2 - 2*x14 + 4*x15 - 2*x16, 
                    2*x4*y3^2 - 4*x5*y3^2 + 2*x6*y3^2 - 2*x6 + 4*x7 - 2*x8}:
        elif test_case ="small_sys_low_deg" then 
            Sys:={x1+y1*x2+y2*x3-1,y2*x1+x2+y1*x3-2,(y1-y2)*x1-x2+y2*x3-7}:
        elif test_case ="small_Sys" then
            Sys:={x1+y1*x2+y1-3,y2*x1+x2+y1-1}:
        elif test_case ="mike" then
            Sys:={y1*x1+y1*x2-1, y1*y2*x1-x2-1}:
        elif test_case = "example" then
            # Your example: 
            Sys := {(y1*y2-1)*x1 + (y1^2-2*y1+3), (y1*y2-1)*x2 + (y1*y2-y1-3*y2+1)}:

        elif test_case = "bsbug" then
            Sys := { (2*y3^2*y4-y3*y4^2+3*y3*y4-y4^2+y3+y4+1)*x1 = y3*y4^2 };

        elif test_case = 1 then
            print("In test_case = 1 ");
            # ff :=x+5*y:
            # gg := x*y+3: Vars := [x, y]:
            # gg:= expand((x-(1/2))*(y-(1/3))*(z-(1/5))): Vars:=[x,y,z]:
            # ff:=x+5*y:
            # gg:=x*y+3:
            ff:=y;
            gg:=x-4;
            Vars:=[x,y]:
            return Vars,ff,gg,numelems(Vars),1,Vars: 
        end if:
    elif nargs > 1 then
        print("In nargs > 1 ");
        print("test_case =", args[1]);
        print("num_var =", args[2]);
        print("num_terms =", args[3]);
        print("den_terms =", args[4]);
        # e.g. data_generator("rand", 5, 7, 9)
        Vars := [ seq( x||i, i = 1..args[2] ) ]:
        if test_case = "rand" then
            ff := randpoly(Vars, terms = args[3]) :
            gg := randpoly(Vars, terms = args[4]) :
            return Vars,ff,gg,numelems(Vars),1,Vars:
        end if;
        if test_case ="rat_rand" then 
            ff:= randpoly(Vars, coeffs=RandRational(args[5]), terms=args[3]);
            gg:= randpoly(Vars, coeffs=RandRational(args[6]), terms=args[4]);


            return Vars,ff,gg,numelems(Vars),1,Vars:    
        end if:
    end if:



    Vars := { seq( x||i, i=1..nops(Sys) )}:
    params := indets(Sys) minus Vars:
    return Sys,Vars,convert(params,list),nops(params),nops(Vars):
end proc:


RandRational := proc(N::posint)
    # Returns a 0-arg proc producing random rationals with bounded numerator/denominator
    return proc() local a,b;
        a := rand(-N..N)();      # numerator (can be zero)
        b := rand(1..N)();       # denominator (never zero)
        if a = 0 then 0 else a/b fi
    end proc;
end proc:
