RCR:=proc(B,num_vars::integer,vars::list,P::list)
    print("In RCR"):
    local Numerator1,Denominator1,Numerator2,Denominator2,Num_coeff_images,Den_coeff_images,temp_n,
    temp_d,M,i,j,num_coeff1,num_mono1,den_coeff1,den_mono1,num_coeff2,num_mono2,den_coeff2,den_mono2,
    term_wise_num_coeff_images,term_wise_den_coeff_images,U_N,U_D,RCN,RCD,num,den:
    for i from 1 to nops(P) do
    num_coeff1,num_mono1,den_coeff1,den_mono1:=MRFI(B,num_var,vars,P[1]):
    num_coeff2,num_mono2,den_coeff2,den_mono2:=MRFI(B,num_var,vars,P[2]):
    print("num_coeff1= ",num_coeff1):
    print("num_mono1= ",num_mono1):
    print("den_coeff1= ",den_coeff1):
    print("den_mono1= ",den_mono1):
    print("num_coeff2= ",num_coeff2):
    print("num_mono2= ",num_mono2):
    print("den_coeff2= ",den_coeff2):
    print("den_mono2= ",den_mono2):
    # Num_coeff_images:=[[coeffs(Numerator1,vars,'mterms')],[coeffs(Numerator2,vars,'mterms')]]:
    # Den_coeff_images:=[[coeffs(Denominator1,vars,'mterms')],[coeffs(Denominator2,vars,'mterms')]]:
    # print("Num_coeff_images= ",Num_coeff_images):
    # print("Den_coeff_images= ",Den_coeff_images):
    term_wise_num_coeff_images:=term_wise_images([num_coeff1,num_coeff2]):
    term_wise_den_coeff_images:=term_wise_images([den_coeff1,den_coeff2]):
    # term_wise_den_coeff_images:=term_wise_images(Den_coeff_images):
    print("term_wise_num_images= ",term_wise_num_coeff_images):
    print("term_wise_den_images= ",term_wise_den_coeff_images):
    U_N:=[]:
    U_D:=[]:
    for i from 1 to nops(term_wise_num_coeff_images) do 
        U_N:=[op(U_N),chrem(term_wise_num_coeff_images[i],P)]:
    end do:
    for i from 1 to nops(term_wise_den_coeff_images) do 
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
    num:= construct_final_polynomial(RCN,num_mono1):
    # print("num= ",num):
    den:= construct_final_polynomial(RCD,den_mono1):
    # print("den= ",den):
    return num,den:
end proc: