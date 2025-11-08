with(IntegerRelations):
RVR:=proc(L,p)
    local M,E,H,B,i,rvr,num_coeff:
    print("In RVR"):
    M:=p:
    print("M = ",M):
    E:=ceil(sqrt(p)):
    print("E = ",E):
    num_coeff:=numelems(L):
    H:=L:
    print("H= ",H):
    H:=[op(H),E]:
    print("H= ",H):
    B := Matrix(num_coeff, num_coeff, shape = identity):
    B:=p*B:
    # print("B before adding zero vector and H= ",B):
    B:=<B|Vector([seq(0,i=1..num_coeff)])>:
    # print("B after adding zero vector= ",B):
    B:=<B;H>:
    print("B after adding H= ",B):
    with(IntegerRelations);
    rvr:=LLL(B):
    print("LLL = ",rvr):
    # v:=rvr[1,1..-2]:
    # q := v[num_coeff+1] / E;
    # a := Vector(num_coeff, i-> v[i]);
    # num_coeff = nops(L) = 2 in your example
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
    end if;

    # verify congruences a_i ≡ q * L[i] (mod M)
    ok := true;
    for i to num_coeff do
        if irem(a[i] - q*L[i], M) <> 0 then ok := false; break; end if;
    od;
    if not ok then error "verification failed"; end if;

    # Return: (rational vector, q, integer numerators)
    # return [seq(a[i]/q, i=1..num_coeff)], q, convert(a, list);
    return [seq(a[i]/q, i=1..num_coeff)]:
end proc:






DecodeRow := proc(v::Vector, X::posint, M::posint, c::Vector)
    local k, a, q, g, i, ok;
    k := LinearAlgebra:-Dimension(v)-1;
    print("k = ",k):
    q := v[k+1] / X;  # should be integer
    print("q = ",q):
    a := Vector[k](i->v[i]);
    print("a = ",a);

    # normalize
    g := igcd(seq(a[i], i=1..k), q);
    print("g = ",g);
    q := q/g;  a := Vector[k](i->a[i]/g);
    print("q = ",q);
    if q < 0 then q := -q; a := -a; fi;

    # verify
    ok := true;
    for i to k do if irem(a[i] - q*c[i], M) <> 0 then ok := false; break; fi; od;
    if not ok then error "verification failed"; fi;

    return [seq(a[i]/q, i=1..k)], q, convert(a, list);
end proc;

# Example with your numbers:
v := Vector([7, 15, 1621935]):
X := 46341:  M := 2147483647:
c := Vector([858993459, 1227133513]):
DecodeRow(v, X, M, c);  # returns [[1/5, 3/7], 35, [7,15