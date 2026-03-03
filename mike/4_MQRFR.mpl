########################################################
# 4. MQRFR - Modified QR for Rational Functions
########################################################

MQRFR:=proc(r0::polynom,r1::polynom,t0::integer,t1::integer,p::prime)
    description "Returns the numerator and denominator polynomials of a rational function":

    # 1. Parameters:
    #   -----------
    # r0 : polynomial 
    #   The product polynomial m(x) whose roots are the evaluation points.
    # r1 : polynomial
    #   The interpolation polynomial u(x) constructed from the evaluation data.
    # t0 : integer
    #   Initialized with 0, used in the Extended Euclidean Algorithm.
    # t1 : integer
    #   Initialized with 1, used in the Extended Euclidean Algorithm.
    # p : prime number
    #     The characteristic of the finite field Z_p
    #     All arithmetic operations are performed modulo p
    #
    # 2. Returns:
    #   --------
    # A tuple (f, g, dq, lcg) where:
    # f : polynomial
    #   The numerator polynomial of the rational function normalized by the leading coefficient of g.
    # g : polynomial
    #   The denominator polynomial of the rational function normalized by its leading coefficient of g.
    # dq : integer
    #   The degree of the denominator polynomial g.
    # lcg : integer
    #   The leading coefficient of the denominator polynomial g.
    #
    # 3. Mathematical Description:
    #   ------------------------
    # This procedure implements a Maximal quotient rational reconstruction algorithm to compute the numerator and denominator
    # of a rational function from its interpolation polynomial and the product polynomial of evaluation points.
    # It uses the Extended Euclidean Algorithm to find polynomials f and g such that:
    #     f(x) ≡ u(x) * g(x) (mod m(x))
    # where u(x) is the interpolation polynomial and m(x) is the product polynomial.
    # The algorithm iteratively computes quotients and updates polynomials until the remainder polynomial becomes zero.
    # The degree of the denominator polynomial g is tracked.
    #
    # 4. Algorithm Complexity:
    #   --------------------
    # Time Complexity: O(d^2) where d is the degree of the polynomial m(x).
    # Space Complexity: O(d) for storing intermediate polynomials.

    # 5. References:
    #   ----------
    # [1] Monagan, M. (2004). "An almost optimal algorithm for rational reconstruction."

 
    print("---------------------------------------------------------"):
    print("In MQRFR"):


    local r,t,q,i,f,g,qmax,lcg:
    r[0]:=r0:
    r[1]:=r1:
    t[0]:=t0:
    t[1]:=t1:
    # lprint("r0=",r0):
    # lprint("r1=",r1):
    # lprint("t0=",t0):
    # lprint("t1=",t1):
    f:=r0:
    g:=t1:
    qmax:=0:
    i:=1:
    while r[i] <> 0 do
        q[i]:= Quo(r[i-1],r[i],x,'r[i+1]') mod p:
        if degree(q[i],x)> qmax then 
            qmax:=degree(q[i],x):
            lprint("q[",i,"]=",q[i]," degree q=",qmax):
            f:=r[i]:
            g:=t[i]:
            # lprint("r[",i-1,"]=",r[i-1]):
            # lprint("q[",i,"]=",q[i]):
            lprint("f=",f):
            lprint("g=",g):
        end if:
        t[i+1]:=Expand(t[i-1]-q[i]*t[i])mod p:
        i:=i+1:
        if qmax <=1 or gcd(f,g) <> 1 or g = 0 then 
            FAIL:
        end if:
    end do:
    lcg:=lcoeff(g):
    # lprint("lcg=",lcg):
    return f/lcg mod p,g/lcg mod p,qmax,lcg :
end proc:
