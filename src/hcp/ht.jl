#=
= Step height move.
=#

#= Structs ====================================================================#

function ht_init(Y, S::State, X::Auxiliary) # hts.
    X.hts[:] = S.hts
end

function ht_unconstrained(X::Auxiliary) # lhts.
    @. X.lhts = log(X.hts)
end

function ht_constrained(X::Auxiliary) # hts.
    @. X.hts = exp(X.lhts)
end

function ht_momenta(C::Constants, X::Auxiliary) # plhts
    X.plhts[:] = randn(C.k+1)
end

#= Move =======================================================================#

function ht_move(Y, C::Constants, S::State, X::Auxiliary, U::Counts)

    ε = .8C.lh + rand()*.4C.lh # base step size ± 20%.
    ht_init(Y,S,X)
    ht_unconstrained(X)
    ht_momenta(C,X)
    ep = ht_potential(C,X) # previous energy; also computes the gradient.
    kp = ht_kinetic(X)

    # ht_println(X) # debug.
    # println("ep: $ep, kp: $kp") # debug.

    ht_pstep(.5ε,X)
    for i = 1:C.Lh
        ht_qstep(ε,X)
        if i != C.Lh
            ht_potential(C,X)
            ht_pstep(ε,X)
        end
    end
    en = ht_potential(C,X)
    ht_pstep(.5ε,X)
    kn = ht_kinetic(X)

    # ht_println(X) # debug.
    # println("en: $en, kn: $kn") # debug.
    # println("diff: $(ep+kp-en-kn)") # debug.
    
    p = rand()
    l = exp(ep+kp-en-kn)
    # println((p < l) ? "$p < $l" : "$p ≥ $l") # debug.

    if p < l
        S.hts[:] = X.hts
        U.haccepts += 1
    else
        ht_init(Y,S,X)      # leave the auxiliary parameters
        ht_unconstrained(X) # in a consistent state.
    end
    U.hprops += 1

end

function ht_pstep(s, X::Auxiliary)
    @. X.plhts -= s * X.∇lhts
end

function ht_qstep(s, X::Auxiliary)
    @. X.lhts += s * X.plhts
    ht_constrained(X)
end

#= Energy =====================================================================#

function ht_potential(C::Constants, X::Auxiliary)
    k,R,α,β = C.k, C.R, C.α, C.β
    hts, widths, counts = X.hts, X.widths, X.counts
    locs = @view(X.locs[2:end-1])

    e = sum(@. (widths+β)*hts - (counts+α)*X.lhts - log(widths)) -
        sum(@. log(locs*(R-locs)))
    @. X.∇lhts = (widths+β)*hts - (counts+α)
    return e
end

ht_kinetic(X::Auxiliary) = .5*sum(X.plhts.^2)

#= I/O ========================================================================#

function ht_println(X::Auxiliary)
    println("X: locs = $(string(X.locs)), hts = $(string(X.hts))\n",
            "   lhts = $(string(X.lhts)), ∇lhts = $(string(X.∇lhts)),",
              " plhts = $(string(X.plhts))")
end

function string(A::Vector{Float64})
    io = IOBuffer()
    write(io, "[")
    for a in @view(A[1:end-1])
        write(io, @sprintf("%g, ", a))
    end
    write(io, @sprintf("%g]", A[end]))
    return String(take!(io))
end