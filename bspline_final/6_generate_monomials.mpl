########################################################
# 6. GENERATE MONOMIALS
########################################################

generate_monomials:=proc(roots_,num_var,prime_points,vars)
    local m,mm,i,j,counter,M_,rem: 
    M_:=Vector(numelems(roots_),0):
    print("In generate_monomials"):
    print("roots_=",roots_):
    
    for i from 1 to numelems(roots_) do
        if(roots_[i]=0)then 
            print("roots_[",i,"] = 0"):
            return FAIL: 
        end if:
        mm:=roots_[i]:
        m:=1:
        for j from 1 to numelems(prime_points) do
            counter:=0:
            while mm mod prime_points[j] = 0 do
                mm:=iquo(mm,prime_points[j],'rem'):
                counter:=counter+1:
            end do:
            m:=m*vars[j]^counter:
        end do:
        M_[i]:=m:
    end do:
    
    if mm<> 1 then 
        print("Warning: mm=",mm," (should be 1)"):
        return FAIL:
    end if:
    
    return convert(M_,list):
end proc:
