data_generator:=proc(test_case,n,num_terms,den_terms,p)
    local var,T,ff,gg,i,r:
    var:=[seq(x||i,i=1..n)];
    
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
        gg:=gg/lcoeff(gg,order=grlex(seq(var[i],i=1..n))) mod p:
    
    elif test_case = "rat_coeff" then
        ff:=(3/4)*x*y+(7/19)*z:
        gg:=(11/47)*x+(2/31)*y*z:
        gg:=gg/lcoeff(gg,order=grlex(x,y,z)):   
        return [x,y,z],T,ff,gg:   
    
    elif test_case = "rat_rand"  then
        #r:=rand(-100..100): 
        r:=rand(1..200):
        ff:=randpoly(var,coeffs=proc() r()/r() end proc, terms =num_terms):
        gg:=randpoly(var,coeffs=proc() r()/r() end proc, terms =den_terms):
        print("ff= ",ff):
        print("gg= ",gg):
        gg:=gg/lcoeff(gg,order=grlex(seq(var[i],i=1..n))):
    end if:
    return var,T,ff,gg:
end proc: