projection_image_phi:=proc(num_var,alpha,beta_,sigma_,p,T)
    print("In projection_image_phi");
    # print("alpha = ",alpha);
    # print("beta_ = ",beta_);
    # print("sigma_ = ",sigma_);   
    # print("num_var = ",num_var);
    # print("T = ",T);
    local phi_,nv,np,i:

        for np from 1 to T do # projection ring_homomorphism
            phi_[np][1]:=alpha[np]:
            # print("phi_[",np,"][1] = ",phi_[np][1]);
            for nv from 2 to num_var do 
                # print("nv = ",nv);
                phi_[np][nv]:=beta_[nv-1]*alpha[np]-beta_[nv-1]*sigma_[1]+sigma_[nv] mod p:
            end do:
        end do:
    return [seq(convert(phi_[i],list),i=1..T)]:
end proc: