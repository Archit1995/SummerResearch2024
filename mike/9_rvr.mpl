with(IntegerRelations):
RVR:=proc(L,p,field)
    local M,E,H,B,i,rvr,num_coeff,a,q,g,ok,qX:
    # print("In RVR"):
    M:=p:
    # print("M = ",M):
    E:=ceil(sqrt(p)):
    # print("E = ",E):
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
    # print("LLL = ",rvr):
    # return rvr[1,1..-2]:

    if field = "Z" then return rvr[1,1..-2]: 

    elif field = "Q" then  # The returns coefficients unto a scalar constant. The constant is lcoeff(OG_Denom).
        qX := rvr[1, num_coeff+1];            
        print("qX = ",qX):
        a  := Vector(num_coeff, i-> rvr[1, i]);    # first num_coeff entries
        print("a = ",a):
        q  := qX / E;
        print("q = ",q):
        # normalize by gcd and sign
        g := igcd(q, seq(a[i], i=1..num_coeff));
        q := iquo(q, g);
        a := Vector(num_coeff, i-> iquo(a[i], g));
        if q < 0 then
            q := -q;
            a := Vector(num_coeff, i-> -a[i]);
        end if:

        # verify congruences a_i ≡ q * L[i] (mod M)
        ok := true:
        for i to num_coeff do
            if irem(a[i] - q*L[i], M) <> 0 then ok := false; break: end if:
        od;
        if not ok then error "verification failed": end if:

        # Return: (rational vector, q, integer numerators)
        # return [seq(a[i]/q, i=1..num_coeff)], q, convert(a, list);
        return [seq(a[i]/q, i=1..num_coeff)]:
    end if:
end proc:






