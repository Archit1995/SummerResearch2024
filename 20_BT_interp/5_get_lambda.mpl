get_lambda:=proc(num_,den_,T,anchor_points,num_var,p)# Needs correction find a way to keep roots in mem for t 
# and t+1 to check if the number of roots is equal to the number of terms. 
    print("In get_num_terms_lambda_roots");
    local t,flag,Y,Lambda,terms,R:
    t:=T:
    i:=1;
    while true do 
        print("i=",i):
        print("t=",t):
        anchor_point_powers:=generate_powers(t,anchor_points,num_var,p):
        print("anchor_point_powers=",anchor_point_powers):
        Y:=[seq(B(prime_powers[p_p],p),p_p=1..numelems(prime_powers))]:
        print("Y=",Y):
        Lambda := BMEA(Y,p,Z):
        print("Lambda=",Lambda):
        terms[i]:=degree(Lambda,Z):
        print("terms=",terms[i]):
        R := Roots(Lambda) mod p:
        print("R=",R):
        print("num of roots of lambda=",nops(R)):
        # if i=1 then 
        if nops(R)<terms[i] then 
            t:=t*2:
        # else 
        elif terms[i-1]=terms[i] and terms[i] = nops(R) then 
            # print("IN TERMINATION oF GET_NUM_TERMS_LAMBDA_ROOTS");
                return terms[i],Lambda,R,Y:
        end if: 
        # end if:
        i:=i+1:
    end do:
    print("terms=",terms[i]):
end proc: