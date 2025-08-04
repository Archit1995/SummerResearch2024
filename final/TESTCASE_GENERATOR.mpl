data_generator:=proc(test_case,n,num_terms,den_terms)
    local var,p,T,ff,gg,i:
    var:=[seq(x||i,i=1..n)];
    p:=2^31-1:
	T:=4:
    if (test_case = 1) then 
        ff:=112*x*y+3*z:
        gg:=x+31*y*z:
    elif(test_case = "denominator_zero") then 
        ff:=x^2*y^2+3:
        # gg:=expand((x-2)*(y-7)*(z-5)) mod p:
        # gg:=expand((x-8)*(x-7)*(y-9)*(y-13)*(z-125)*(z-121)*(z-17)) mod p:# Works!!
        # gg:=expand((x-8)*(y-27)*(z-125)) mod p:# Works!!
        gg:=expand((x-2)*(y-3)*(z-5)) mod p:
    elif test_case = "numerator_zero" then 
        ff:=expand((x-2)*(y-3)*(z-5));
        gg:=x^2+y^2+3:    
    elif test_case = "rand"  then
        ff:=randpoly(var,terms =num_terms) mod p;
        gg:=randpoly(var,terms =den_terms) mod p;
        # gg:=gg/lcoeff(gg,order=grlex(seq(var[i],i=1..n))) mod p:
    
    elif (test_case = 99) then
       ff:=randpoly([x,y,z],degree =99) mod p;
        gg:=randpoly([x,y,z],terms =7) mod p;
        gg:=gg/lcoeff(gg,order=grlex(x,y)) mod p:      
    end if:
    return var,p,T,ff,gg:
end proc: