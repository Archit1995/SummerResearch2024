########################################################
# 5. BMEA - Berlekamp-Massey Extended Algorithm
########################################################

BMEA := proc(v::list,p::posint,Z::name)
    description "Implements the Berlekamp-Massey Extended Algorithm to find the minimal polynomial of a given sequence":
# 1. Parameters:
# -----------
# v : list
#     A list of integers representing the sequence for which the minimal polynomial is to be found.
# p : positive integer
#     A prime number representing the modulus for arithmetic operations.
# Z : name
#     A variable name used in polynomial expressions.
# 2. Returns:
# --------
# polynom
#     The minimal polynomial of the input sequence v over the finite field Z_p.
# 3. Mathematical Description:
# ------------------------
# The Berlekamp-Massey Extended Algorithm is an efficient method for finding the minimal polynomial
# of a given sequence. The minimal polynomial is the polynomial of least degree that generates the sequence
# when evaluated at successive integers. The algorithm iteratively updates the candidate polynomial
# based on discrepancies observed in the sequence, ensuring that the polynomial remains minimal.
# 4. Algorithm Complexity:
# --------------------
# Time Complexity: O(n^2) where n is the length of the input sequence v.
# Space Complexity: O(n) for storing intermediate polynomials.  
# 5. References:
# ----------
# [1] Berlekamp, E. R. (1968). "Algebraic Coding Theory."
# [2] Massey, J. L. (1969). "Shift-register synthesis and BCH decoding."


    local n,m,R0,R1,V0,V1,i,Q:
    print("______________________________________________________________________________");
    print("In BMEA");
    lprint("v=",v);    
    n := iquo( nops(v), 2 ):
    # lprint("n=",n);
    m := 2*n-1:
    # lprint("m=",m);
    R0 := Z^(2*n):
    # lprint("R0=",R0);
    R1 := add( v[m+1-i]*Z^i, i=0..m ) mod p:
    # lprint("R1=",R1);
    V0 := 0:
    V1 := 1:
    while n <= degree(R1,Z) do
        R0,R1 := R1,Rem(R0,R1,Z,'Q') mod p:
        # lprint("R0=",R0);   
        # lprint("R1=",R1);

        V0,V1 := V1,Expand(V0-Q*V1) mod p:
        # lprint("V0=",V0);   
        # lprint("V1=",V1);
    od:
    i := 1/lcoeff(V1,Z) mod p:
    return i*V1 mod p:
end:
