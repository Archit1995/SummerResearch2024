term_wise_images:=proc(coeff_images::list)
    print("In term_wise_images"):
    print("coeff_images= ",coeff_images):
    local i,temp:
    temp:=[]:
            for i from 1 to nops(coeff_images[1]) do 
                temp:=[op(temp),[coeff_images[1][i],coeff_images[2][i]]]:
        end do:
    print("temp= ",temp):
    return temp:
end proc: