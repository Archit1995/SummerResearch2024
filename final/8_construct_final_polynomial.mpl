construct_final_polynomial:=proc(coeff_,Monomials)
    local i,f,n:
    f:=0:
    for i from 1 to numelems(coeff_) do
        f:=f+coeff_[i]*Monomials[i]:
    end do:
    return f:
end proc: