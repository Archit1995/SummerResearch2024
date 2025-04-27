data_generator:=proc(test_case,deg,p)
    local num_var,vars,T,ff,gg,r:
    num_deg:=deg[1];
    den_deg:=deg[2];
    r:=rand(p):
    T:=4:
	# p:=811:
	# p:=19:
    ff:=randpoly(x,degree=deg[1]) mod p:
    gg:=randpoly(x,degree=deg[2]) mod p:
    if (test_case = 1) then #working 
        actual_lc:=lcoeff(gg):
    elif (test_case = 2) then 
        gg:=gg*r() mod p:
        actual_lc:=lcoeff(gg):
    else 
        actual_lc:=lcoeff(gg):
        ff:=ff/actual_lc mod p:
        gg:=gg/actual_lc mod p:
    end if:
    return ff,gg,actual_lc:
end proc: