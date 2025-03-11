# with(Bits):
# generate_monomials:=proc(roots_,num_var,prime_points,vars)# needs correction- We are getting roots of lambda polynomial that are not 
#     # factors of 2,3 and 5;
#     print("In generate_monomials");
#     print("roots_=",roots_);
#     local m,mm,i,j,counter,M_: 
#     M_:=Vector(numelems(roots_),0):
#     print("r=",numelems(roots_)):
#     print("vars=",vars):  
#     print("vars[1]=",vars[1]):   
#     print("prime_points=",prime_points):     
#     check:=roots_[1][1]:
#     for ii from 2 to numelems(roots_[1]) do
#         check:=Bits[And](check,roots_[1][ii]):
#         end do;
#     print("check=",check):
#     if check=1 then start:= 2 else start:=1 end if;
#     for i from start to  numelems(roots_) do # each root
#         print("i=",i):

#         m:=1:
#         print("roots_[i]=",roots_[i]):
#         for j from 1 to numelems(prime_points) do #  each prime
#             print("j=",j):
#             counter:=0:
#             print("j=",j):
#             mm:=roots_[i][j]:
#             print("mm=",mm):
#             while mm mod prime_points[j] = 0 do #repeated division
#                 print("in while loop"); 
#                 print("prime_points[j]=",prime_points[j]):
#                 mm:=iquo(mm,prime_points[j]):
#                 print("mm=",mm):
#                 counter:=counter+1:
#                 print("counter=",counter):
#                 print("================================================"):
#             end do:
#             m:=m*vars[j]^counter:# each monomial
#             print("m=",m):
#             print("-----------------------------------------------------------------"):
#         end do:
#         print("m=",m):
#         print("i=",i):
#         M_[i]:=m:
#         print("M[i]=",M_[i]):
#         print("______________________________________"):
#     end do:
#     print([seq(ifactor(roots_[i]),i=1..numelems(roots_))]);# We are getting roots of lambda polynomial that are not 
#     # factors of 2,3 and 5;
#     return convert(M_,list):
# end proc:
generate_monomials:=proc(roots_,num_var,prime_points,vars)# needs correction- We are getting roots of lambda polynomial that are not 
    # factors of 2,3 and 5;
    local m,mm,i,j,counter,M_: 
    M_:=Vector(numelems(roots_),0):
    print("r=",numelems(roots_)):
    for i from 1 to  numelems(roots_) do # each root
        print("i=",i):
        mm:=roots_[i]:
        m:=1:
        print("roots_[i]=",roots_[i]):
        for j from 1 to numelems(prime_points) do #  each prime
            counter:=0:
            print("j=",j):
            while mm mod prime_points[j] = 0 do #repeated division
                print("prime_points[j]=",prime_points[j]):
                mm:=iquo(mm,prime_points[j]):
                print("mm=",mm):
                counter:=counter+1:
                print("counter=",counter):
                print("================================================"):
            end do:
            m:=m*vars[j]^counter:# each monomial
            print("m=",m):
            print("-----------------------------------------------------------------"):
        end do:
        print("m=",m):
        print("i=",i):
        M_[i]:=m:
        print("M[i]=",M_[i]):
        print("______________________________________"):
    end do:
    print([seq(ifactor(roots_[i]),i=1..numelems(roots_))]);# We are getting roots of lambda polynomial that are not 
    # factors of 2,3 and 5;
    return convert(M_,list):
end proc: