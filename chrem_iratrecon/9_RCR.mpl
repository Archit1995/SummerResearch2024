RCR:=proc(B,num_vars::integer,vars::list,P::list)
    print("In RCR"):
    local Numerator1,Denominator1,Numerator2,Denominator2,Num_coeff_images,Den_coeff_images,temp_n,
    temp_d,M,i,j:
    Numerator1,Denominator1:=MRFI(B,num_var,vars,P[1]):
    Numerator2,Denominator2:=MRFI(B,num_var,vars,P[2]):
    Num_coeff_images:=[[coeffs(Numerator1,vars,'mterms')],[coeffs(Numerator2,vars,'mterms')]]:
    Den_coeff_images:=[[coeffs(Denominator1,vars,'mterms')],[coeffs(Denominator2,vars,'mterms')]]:
    print("Num_coeff_images= ",Num_coeff_images):
    print("Den_coeff_images= ",Den_coeff_images):
    term_wise_num_coeff_images:=term_wise_images(Num_coeff_images):
    term_wise_den_coeff_images:=term_wise_images(Den_coeff_images):
    print("term_wise_num_images= ",term_wise_num_coeff_images):
    print("term_wise_den_images= ",term_wise_den_coeff_images):
    U_N:=[]:
    U_D:=[]:
    for i from 1 to nops(term_wise_num_coeff_images) do 
        U_N:=[op(U_N),chrem(term_wise_num_coeff_images[i],P)]:
        U_D:=[op(U_D),chrem(term_wise_den_coeff_images[i],P)]:
    end do:
    print("U_N= ",U_N):
    print("U_D= ",U_D):
    print("P= ",P):
    M:=product(P[j],j=1..nops(P)):
    print("M= ",M):
    RCN:=[seq(iratrecon(U_N[i],M),i=1..nops(U_N))]:
    RCD:=[seq(iratrecon(U_D[i],M),i=1..nops(U_D))]:
    print("RCN= ",RCN):
    print("RCD= ",RCD):
end proc: