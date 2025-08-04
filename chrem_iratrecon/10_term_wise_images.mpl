term_wise_images:=proc(coeff_images::list)
    local i,j,temp,temp_termwise:
    temp_termwise:=[]:
        for i from 1 to nops(coeff_images) do 
            temp:=[]:
            for j from 1 to nops(coeff_images[1]) do 
                temp:=[op(temp),coeff_images[j][i]]:
            end do:
            temp_termwise:=[op(temp_termwise),temp]:
        end do:
    return temp_termwise:
end proc: