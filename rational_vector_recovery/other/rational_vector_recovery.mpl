# Rational Vector Recovery Algorithm from Kaltofen-Yang 2007
# Based on Algorithm on page 208

RationalVectorRecovery := proc(H_list::list, M::posint, s::posint:=5, max_trials::posint:=10)
    local t, trial, gamma, H, cf_list, l, U, V, E, delta, G_candidates, G, D, i, Fi, all_good;
    
    description "Recover rational vector from modular images using continued fractions";
    
    t := nops(H_list):
    
    # Try multiple random linear combinations
    for trial from 1 to max_trials do
        print("Trial", trial):
        
        # Step vr1: Random linear combination
        gamma := [seq(rand(-s..s)(), i=1..t)]:
        H := add(gamma[i] * H_list[i], i=1..t) mod M:
        
        if H = 0 then
            print("  H=0, trying next trial"):
            next:
        end if:
        
        print("  Random combination H =", H):
        
        # Step vr2: Compute continued fraction approximations of H/M
        cf_list := cfrac(H, M, 'quotients'):
        
        # Try each convergent
        for l from 1 to nops(cf_list) do
            # Get the l-th convergent U/V
            V := denom(cf_list[l]):
            U := numer(cf_list[l]):
            
            # Step vr3: Set E and compute bounds
            E := V + 1:
            
            # Check GCD(H, M) < E
            delta := igcd(H, M):
            if delta >= E then
                next:  # Go to next convergent
            end if:
            
            # Compute maximum bound D satisfying (D-1)(E-1) < M < D*E
            # Solving for D: D > M/E and D <= (M + E - 1)/E
            D := iquo(M + E - 1, E):
            
            # Compute G candidates using Theorem 4.1
            G_candidates := []:
            
            # Primary solution: G1 = V
            G := V:
            if igcd(G, M) = 1 then
                G_candidates := [op(G_candidates), G]:
            end if:
            
            # There might be a second solution G2 (see Theorem 4.1)
            # For simplicity, we only use G1 = V here
            
            # Step vr4: Test each candidate G
            for G in G_candidates do
                all_good := true:
                
                # Check if |G*Hi mod M| < D for all i
                for i from 1 to t do
                    Fi := mods(G * H_list[i], M):  # Symmetric remainder
                    
                    if abs(Fi) >= D then
                        all_good := false:
                        break:
                    end if:
                end do:
                
                if all_good then
                    print("  SUCCESS! Found G =", G, "D =", D):
                    # Return the common denominator G and bound D
                    return G, D:
                end if:
            end do:
        end do:
    end do:
    
    # If we get here, all trials failed
    print("FAILURE: Could not find valid G after", max_trials, "trials"):
    return FAIL:
end proc:


# Helper function to compute continued fraction convergents
# Maple's cfrac returns the convergents directly
cfrac := proc(n, d, opt)
    local q, r, convs, i, p0, p1, p2, q0, q1, q2, quotients;
    
    # Compute continued fraction quotients
    quotients := []:
    p0 := n:
    q0 := d:
    
    while q0 <> 0 do
        q := iquo(p0, q0, 'r'):
        quotients := [op(quotients), q]:
        p0 := q0:
        q0 := r:
    end do:
    
    # Compute convergents from quotients
    convs := []:
    p0 := 0: q0 := 1:  # P_{-1}/Q_{-1}
    p1 := 1: q1 := 0:  # P_0/Q_0
    
    for i from 1 to nops(quotients) do
        p2 := quotients[i] * p1 + p0:
        q2 := quotients[i] * q1 + q0:
        convs := [op(convs), p2/q2]:
        p0 := p1: q0 := q1:
        p1 := p2: q1 := q2:
    end do:
    
    if nargs > 2 and opt = 'quotients' then
        return convs:
    else
        return convs:
    end if:
end proc:


# Function to recover the actual rational numbers once G and D are known
RecoverRationalVector := proc(H_list::list, M::posint, G::posint, D::posint)
    local i, Fi, rational_vector:
    
    rational_vector := []:
    
    for i from 1 to nops(H_list) do
        Fi := mods(G * H_list[i], M):
        
        if abs(Fi) >= D then
            error "Numerator bound violated for component %1", i:
        end if:
        
        # The rational number is Fi/G
        rational_vector := [op(rational_vector), Fi/G]:
    end do:
    
    return rational_vector:
end proc:


# Example usage
TestRationalVectorRecovery := proc()
    local V, M, H, G, D, recovered, i:
    
    print("=== Test Case from Paper (Example 4.2) ==="):
    
    # The true rational vector
    V := [103/5003, 1847/5003, -339/5003, -3772/5003, 1060/5003, 2234/5003, 3085/5003, 4826/5003]:
    
    # Modulus
    M := 2^17:  # = 131072
    
    # Compute modular images
    H := [seq(numer(V[i]) * (1/denom(V[i]) mod M) mod M, i=1..nops(V))]:
    
    print("True vector V =", V):
    print("Modulus M =", M):
    print("Modular images H =", H):
    print(""):
    
    # Recover using the algorithm
    G, D := RationalVectorRecovery(H, M, 5, 10):
    
    if G <> FAIL then
        print(""):
        print("Common denominator G =", G):
        print("Bound D =", D):
        
        recovered := RecoverRationalVector(H, M, G, D):
        print("Recovered vector =", recovered):
        
        # Verify
        print(""):
        print("Verification:");
        for i from 1 to nops(V) do
            if recovered[i] = V[i] then
                print("  Component", i, ": CORRECT"):
            else
                print("  Component", i, ": WRONG - expected", V[i], "got", recovered[i]):
            end if:
        end do:
    end if:
end proc:

TestRationalVectorRecovery();