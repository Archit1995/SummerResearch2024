generate_monomials:=proc(roots_,num_var,prime_points,vars)# needs correction- We are getting roots of lambda polynomial that are not 
    # factors of 2,3 and 5;
    local m,mm,i,j,counter,M_,rem: 
    M_:=Vector(numelems(roots_),0):
    print("In generate_monomials"):
    print("roots_=",roots_):
    print("vars=",vars):
    print("prime_points=",prime_points):
    print("r=",numelems(roots_)):
    for i from 1 to  numelems(roots_) do # each root
        print("i=",i):
        mm:=roots_[i]:
        m:=1:
        print("roots_[i]=",roots_[i]):
        for j from 1 to numelems(prime_points) do #  each prime
            counter:=0:
            print("prime=",j):
            while mm mod prime_points[j] = 0 do #repeated division
                print("sigma=",prime_points[j]):
                mm:=iquo(mm,prime_points[j],'rem'):
                print("quotient =",mm):
                print("remainder =",rem):
                counter:=counter+1:
                print("counter=",counter):
                print("================================================"):
            end do:
            m:=m*vars[j]^counter:# each monomial
            print("m=",m):
            print("-----------------------------------------------------------------"):
        end do:
        # print("m=",m):
        # print("i=",i):
        M_[i]:=m:
        print("M[i]=",M_[i]):
        print("______________________________________"):
    end do:
    print("mm=",mm):#Incase Monomial(2^i,3^i,5,^i) > p, the roots will not be a factor of 2,3 and 5;
    if mm<> 1 then 
        return FAIL:
    end if:
    # For cross checking
    print([seq(ifactor(roots_[i]),i=1..numelems(roots_))]);# We are getting roots of lambda polynomial that are not 
    # factors of 2,3 and 5;
    return convert(M_,list):
end proc: