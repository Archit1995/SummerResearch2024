data_generator:=proc(test_case)
    local num_var,vars,p,T,ff,gg:
    num_var:=3:
    vars:={x,y,z}:
	p:=2^31-1:
	# p:=19:
	T:=4:
	# p:=811:
	# p:=19:
    if (test_case = 1) then #working 
        ff:=112*x*y+3*z:
        gg:=x+31*y*z:
    elif (test_case = 2) then #working for ff:= x^2*y^2+3, we do not get a monic denominator- Why?
        ff:=x^2*y^2+3:
        gg:=x^2+y^2+3:
    # Degree of num = 2x deg of den
    elif (test_case = 3) then #NOT working num_evals are not matching -Why? The lambda polynomial has no roots
    # The ratio of raw_num to true_num = raw_den to true_den but != f(0)
    # Also, f(0) is different for each evaluation "
    #  RATIO DEPENDS OF SHIFT BETA and NOTon sigma =[2,3,5,7,...]
        ff:=x^2*y^2+x+y+1:# 4 terms 
        gg:=x^2+y^2+x+y+3:# 5 terms
    elif (test_case = 4) then
        ff:=x^3*y^3+x^2*y^2+x*y+1:
        gg:=x^2+y^2+x+y+3:
    elif (test_case = 5) then
        ff:=12*x^3*y^3+15415*x^2*y^2+x*y+87548:
        gg:=x^2+y^2+4586*x*y+89456*x+12413*y+9367823:
    elif (test_case = 6) then # works for small prime [2,3] and not for [7,11]
        ff:=x^6*y^7+x^2*y^2+x*y+1:
        gg:=x^3+y^4+x^2+y+3:
    elif (test_case = 7) then
        ff:=x^12*y^12+x+y+1:# 4 terms 
        gg:=x^6+y^6+x+y+3:# 5 terms
    elif (test_case = 8) then# incorrect answer
        ff:=x^23*y^13+x^6*y^7+x*y+1:
        gg:=x^12+y^32+x^4+y^7+3:            
    elif (test_case = 9) then
        ff:=randpoly([x,y,z],terms =31) mod p;
        gg:=randpoly([x,y,z],terms =7) mod p;
        gg:=gg/lcoeff(gg,order=grlex(x,y)) mod p:
    elif (test_case = 10) then
        ff:=x^2+y^2+89456*x+12413*y+9367823:
        gg:=12*x^3*y^3+15415*x^2*y^2+x*y+87548
    elif (test_case = 44) then
        ff:=x^2+y^2+x+y+3:
        gg:=x^3*y^3+x^2*y^2+x*y+1:
    elif (test_case = 11) then
        ff:=x+y:
        gg:=x*y+1:           
    elif (test_case = 13) then
        ff:=x^3*y^3+x^2*y^2+x*y+1:
        gg:=x^2+y^2+x+y+3:
    elif (test_case = 14) then
        ff:=x^3*y^3+x^2*y^2+x*y+1:
        gg:=x^2+y^2+x+y+3:        
    else
        ff:=x*y+1:
        gg:=x+y:
    end if:
    return num_var,vars,p,T,ff,gg:
end proc: