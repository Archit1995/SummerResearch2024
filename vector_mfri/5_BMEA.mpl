########################################################
# 5. BMEA - Berlekamp-Massey Extended Algorithm
########################################################

BMEA := proc(v::list,p::posint,Z::name)
    local n,m,R0,R1,V0,V1,i,Q:
    print("______________________________________________________________________________");
    print("In BMEA");
    lprint("v=",v);    
    n := iquo( nops(v), 2 ):
    lprint("n=",n);
    m := 2*n-1:
    lprint("m=",m);
    R0 := Z^(2*n):
    lprint("R0=",R0);
    R1 := add( v[m+1-i]*Z^i, i=0..m ) mod p:
    lprint("R1=",R1);
    V0 := 0:
    V1 := 1:
    while n <= degree(R1,Z) do
        R0,R1 := R1,Rem(R0,R1,Z,'Q') mod p:
        lprint("R0=",R0);   
        lprint("R1=",R1);

        V0,V1 := V1,Expand(V0-Q*V1) mod p:
        lprint("V0=",V0);   
        lprint("V1=",V1);
    od:
    i := 1/lcoeff(V1,Z) mod p:
    return i*V1 mod p:
end:
