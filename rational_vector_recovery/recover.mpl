with(IntegerRelations):
RVR:=proc(L,p)
    print("In RVR"):
    M:=p:
    E:=ceil(sqrt(p)):
    num_coeff:=numelems(L):
    H:=L:
    # print("H= ",H):
    H:=[op(H),E]:
    # print("H= ",H):
    B := Matrix(num_coeff, num_coeff, shape = identity):
    B:=p*B:
    # print("B before adding zero vector and H= ",B):
    B:=<B|Vector([seq(0,i=1..num_coeff)])>:
    # print("B after adding zero vector= ",B):
    B:=<B;H>:
    # print("B after adding H= ",B):
    with(IntegerRelations);
    rvr:=LLL(B):
    return rvr[1,1..-2];
end proc: