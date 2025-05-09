data_generator:=proc(test_case)
    local num_var,vars,p,T,ff,gg:
    num_var:=2:
    vars:={x,y}:
	p:=2^31-1:
	# p:=19:
	T:=4:
	# p:=811:
	# p:=19:
    if (test_case = 1) then #working 
        ff:=x*y+1:
        gg:=x+y:
    elif (test_case = 2) then #working DOESNT WORK for ff:= x^2*y^2+3
        ff:=x^2*y^2+1:
        gg:=x^2+y^2+3:
    elif (test_case = 3) then #NOT working
        ff:=x^2*y^2+x+y+1:
        gg:=x^2+y^2+x+y+3:
    elif (test_case = 4) then
        ff:=x^3*y^3+x^2*y^2+x*y+1:
        gg:=x^2+y^2+x+y+3:
    elif (test_case = 5) then
        ff:=x^3*y^3+x^2*y^2+x*y+1:
        gg:=x^2+y^2+x+y+3:
    else
        ff:=x*y+1:
        gg:=x+y:
    end if:
    return num_var,vars,p,T,ff,gg:
end proc: