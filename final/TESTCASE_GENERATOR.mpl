data_generator := proc(test_case)
    local var, p, T, ff, gg, i;
    print("In data_generator");
    print("nargs =", nargs);
    print("test_case =", test_case);
    p := 2^31 - 1:
    T := 4:
    if nargs = 1 then
        var := [x, y, z]:
        print("In nargs = 1 ");
        if test_case = "bspline_small_sys_low_deg1" then 
            ff:=-2*y1*y2^3+y2^3+7*y1^2-9*y2^2+y1:
            gg:=-y1*y2^6-y2^5+y2^3+y1^2-y2^2+y1:
            var:=[y1,y2]:
        elif test_case = "bspline_small_sys_low_deg2" then
            ff:=y2^6-7*y2^2-2*y2^3+2*y2^2+6*y1:
            gg:=-y1*y2^6-y2^5+y2^3+y1^2-y2^2+y1:
            var:=[y1,y2]:
        elif test_case = "bspline_small_sys_low_deg3" then
            ff:=-7*y1*y2^3-y2^3-2*y1-8:
            gg:=-y1*y2^6-y2^5+y2^3+y1^2-y2^2+y1:
            var:=[y1,y2]:
        elif test_case = 1 then
            print("In test_case = 1 ");
            ff := x*y + 3*z:
            gg := x + y*z:
        elif test_case = "denominator_zero" then
            ff := x^2*y^2 + 3:
            gg := expand((x-2)*(y-3)*(z-5)) mod p:
        elif test_case = "numerator_zero" then
            ff := expand((x-2)*(y-3)*(z-5)):
            gg := x^2 + y^2 + 3:
        elif test_case = 99 then
            ff := randpoly([x, y, z], degree = 99) mod p:
            gg := randpoly([x, y, z], terms = 7) mod p:
        end if:
    elif nargs > 1 then
        print("In nargs > 1 ");
        print("test_case =", args[1]);
        print("num_var =", args[2]);
        print("num_terms =", args[3]);
        print("den_terms =", args[4]);
        # e.g. data_generator("rand", 5, 7, 9)
        var := [ seq( x||i, i = 1..args[2] ) ]:
        if test_case = "rand" then
            ff := randpoly(var, terms = args[3]) mod p:
            gg := randpoly(var, terms = args[4]) mod p:
        end if;
    end if:

    # **THIS** returns a 5-value sequence, so you can unpack it:
    return( var, p, T, ff, gg ):
end proc:
