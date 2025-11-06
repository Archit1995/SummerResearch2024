# ------------------------------------------------------------
# Rational vector reconstruction with a common denominator
# from multiple modular image vectors (pairwise coprime moduli).
#
# Input:
#   Hmods :: list(list(integer))  -- list of residue vectors; each element is a list [h1,...,ht] modulo M_k
#   Ms    :: list(integer)         -- list of pairwise-coprime moduli [M1,...,Mm], same length as Hmods
# Options:
#   s        :: nonneg int (default 5)     -- random coeff range: -s..s (vr1)
#   trials   :: pos int (default 40)       -- max random trials (outer loop)
#   seed     :: integer or none            -- optional: to make randomness reproducible
#
# Output:
#   record with fields:
#     num    :: Vector(integer)  -- integer numerators F_1..F_t
#     den    :: integer          -- positive common denominator G
#     rat    :: Vector(rational) -- simplified rational vector F_i/G
#     M      :: integer          -- combined modulus
# or raises error if reconstruction fails.
# ------------------------------------------------------------

RationalVectorRecover := proc(Hmods::list, Ms::list; s::nonnegint:=5, trials::posint:=40, seed::anything:=NULL)
    option autoload;
    uses LinearAlgebra, numtheory;

    local m, t, k, i, CRTVec, symres, M, H, ok, rndCoeffs, Hcomb, CF, convs, l,
          E, D1, D2, D, Gcands, G, gcdGM, testOK, Fi, Fvec, g, gAll, lz, den, numV, ratV,
          rseed, Mi, Hi_list, H_i, Hk, crt_one, gcd_all, simplify_sign, reduce_common;

    # ---------- helpers ----------

    # Symmetric residue in (-M/2, M/2]
    symres := proc(a::integer, M::posint) local r;
        r := irem(a, M);
        if r > floor(M/2) then r := r - M fi;
        return r
    end proc;

    # Chinese remainder for a vector across pairwise coprime moduli
    CRTVec := proc(Hmods::list, Ms::list)
        local m, t, i, k, H, M, curr, crt_pair;
        m := nops(Ms);
        if m <> nops(Hmods) then error "Length mismatch: Hmods vs Ms" fi;
        t := nops(Hmods[1]);
        for k to m do
            if nops(Hmods[k]) <> t then
                error "All residue vectors must have the same length"
            fi
        od;
        # combine componentwise using chrem
        H := Vector(t);
        M := 1;
        for i to t do
            curr := [Hmods[1][i], Ms[1]];
            for k from 2 to m do
                curr := [chrem(curr[1], curr[2], Hmods[k][i], Ms[k]), curr[2]*Ms[k]];
            od;
            H[i] := curr[1];
            M := curr[2];  # final product
        od;
        return H, M
    end proc;

    # Build random small coefficients in [-s..s], not all zero
    rndCoeffs := proc(t::posint, s::nonnegint)
        local gam, j, allzero;
        gam := Vector(t);
        allzero := true;
        for j to t do
            gam[j] := rand(-s..s)();
            if gam[j] <> 0 then allzero := false fi;
        od;
        if allzero then
            # force one random nonzero
            gam[1] := max(1, s);  # simple nudge
        fi;
        return gam
    end proc;

    # Produce the convergents (U_l / V_l) of H/M via continued fractions
    # Returns list of pairs [U,V] with V > 0
    convs := proc(H::integer, M::posint)
        local cf, CV, j, U, V, L;
        # Maple gives a finite/infinite CF; we take convergents until stabilized
        cf := convert(H/M, confrac);
        CV := numtheory[convergents](cf);
        L := [];
        for j to nops(CV) do
            U := numer(CV[j]); V := denom(CV[j]);
            if V > 0 then L := [op(L), [U,V]] fi;
        od;
        return L
    end proc;

    # Optionally reduce common gcd across all numerators and the denominator
    reduce_common := proc(F::Vector, G::posint)
        local gAll, i;
        gAll := G;
        for i to LinearAlgebra[Dimension](F) do
            gAll := igcd(gAll, abs(F[i]));
        od;
        if gAll > 1 then
            for i to LinearAlgebra[Dimension](F) do F[i] := iquo(F[i], gAll) od;
            G := iquo(G, gAll);
        fi;
        # enforce positive denominator
        if G < 0 then
            G := -G;
            for i to LinearAlgebra[Dimension](F) do F[i] := -F[i] od;
        fi;
        return F, G
    end proc;

    # ---------- validate input ----------
    m := nops(Ms);
    if m = 0 then error "Empty moduli list" fi;
    t := nops(Hmods[1]);
    for k to m do
        if type(Ms[k], posint) = false then error "Moduli must be positive integers" fi;
        if nops(Hmods[k]) <> t then error "All residue vectors must have same length" fi;
        # sanity: residues mod Mk
        for i to t do
            if type(Hmods[k][i], integer) = false then error "Residues must be integers" fi;
        od;
    od;

    if seed <> NULL then randomize(seed) fi;

    # ---------- combine all moduli with CRT (componentwise) ----------
    H, M := CRTVec(Hmods, Ms);

    # ---------- outer randomized trials (vr1) ----------
    for lz from 1 to trials do

        # vr1: random linear combination Hcomb = sum gamma_i * H_i  (mod M)
        #      (keep t small? we use all entries; this works well in practice)
        local gamma;
        gamma := rndCoeffs(t, s);
        Hcomb := 0;
        for i to t do
            Hcomb := irem(Hcomb + gamma[i]*H[i], M);
        od;
        if Hcomb = 0 then next fi;  # unlucky, retry

        # vr2: generate convergents of Hcomb / M
        local C; C := convs(Hcomb, M);

        # iterate convergents: for each [U_l, V_l], let E := V_l + 1
        for l to nops(C) do
            local Ul, Vl;
            Ul := C[l][1]; Vl := C[l][2];
            E  := Vl + 1;
            if E <= 1 then next fi;

            # vr3: if gcd(Hcomb, M) >= E, skip (trial fails early)
            if igcd(Hcomb, M) >= E then next fi;

            # find the maximal integer D satisfying (D-1)(E-1) < M < D*E.
            # inequalities => D > M/E and D-1 < M/(E-1)
            D1 := iquo(M,   E) + 1;         # ensures M < D*E
            D2 := iquo(M-1, E-1) + 1;       # ensures (D-1)(E-1) < M
            D  := max(D1, D2);

            # Candidate denominators: G1 = Vl; (optionally also G2 from neighbors per Thm 4.1)
            Gcands := [Vl];

            # (Optional) attempt neighbor-based second candidate if helpful:
            # In many practical cases, G1 suffices thanks to vr4 self-correction.
            # To cover the rare two-solution case you can uncomment below:
            # if l>=2 then Gcands := [op(Gcands), C[l-1][2]] fi;
            # if l+1<=nops(C) then Gcands := [op(Gcands), C[l+1][2]] fi;

            # test each candidate G (vr4)
            for G in Gcands do
                gcdGM := igcd(G, M);
                if gcdGM <> 1 then next fi;

                # test bound | G * H_i (mod M) | < D for all components (symmetric residue)
                testOK := true;
                Fvec   := Vector(t);

                for i to t do
                    Fi := symres( irem(G*H[i], M), M );
                    if abs(Fi) >= D then
                        testOK := false;
                        break
                    fi;
                    Fvec[i] := Fi;
                od;

                if testOK then
                    # Success: normalize (divide common gcd across F and G; make G > 0)
                    Fvec, G := reduce_common(Fvec, G);

                    # Construct outputs
                    den  := G;
                    numV := Vector(t, i->Fvec[i]);
                    ratV := Vector(t, i-> numV[i] / den );

                    return Record(
                        num = numV,
                        den = den,
                        rat = ratV,
                        M   = M
                    );
                fi;  # end success
            od; # next G
        od; # next convergent
    od; # next trial

    error "RationalVectorRecover: failed to reconstruct after %d trials; try increasing modulus or trials", trials;
end proc;


# ------------------------------------------------------------
# Convenience wrapper: accepts raw (vectors, modulus) pairs and calls the core
# Example signature:
#   RationalVectorRecoverFromPairs( [[h^1_1,...,h^1_t], M1],
#                                   [[h^2_1,...,h^2_t], M2], ...,
#                                   options ... )
# ------------------------------------------------------------
RationalVectorRecoverFromPairs := proc()
    local L, Hmods, Ms, i, opts, pos, arg;
    L := [args];
    # split positional pairs and options
    Hmods := []; Ms := []; pos := 1;
    while pos <= nops(L) and type(L[pos], list) and nops(L[pos])=2 do
        Hmods := [op(Hmods), L[pos][1]];
        Ms    := [op(Ms),    L[pos][2]];
        pos := pos+1;
    od;
    opts := op(L[pos..-1]);
    return RationalVectorRecover(Hmods, Ms, opts);
end proc;
