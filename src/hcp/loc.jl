#=
= Step location move.
=#

#= Auxiliary ==================================================================#

function loc_init(Y, S::State, X::Auxiliary) # locs, widths, counts.
    X.locs[:] = S.locs
    loc_metrics(Y,X)
end

function loc_metrics(Y, X::Auxiliary) # widths, counts.
    sp = 0.
    y = 1 # index of the first uncounted datum.
    for i = eachindex(X.widths)
        s = X.locs[i+1]
        X.widths[i] = s-sp
        n = searchsortedfirst(@view(Y[y:end]), s) - 1
        X.counts[i] = n
        sp,y = s,y+n
    end
end

function loc_unconstrained(C::Constants, X::Auxiliary) # tlocs, texps, dlocs.
    locs = @view(X.locs[2:end-1])
    @. X.tlocs = log(locs/(C.R-locs))
end

function loc_constrained(Y, C::Constants, X::Auxiliary) # locs, widths, counts.
    @. X.locs[2:end-1] = C.R/(1+exp.(-X.tlocs))
    loc_metrics(Y,X)
end

function loc_momenta(C::Constants, X::Auxiliary) # ptlocs.
    X.ptlocs[:] = randn(C.k)
end

#= Move =======================================================================#

function loc_move(Y, C::Constants, S::State, X::Auxiliary, U::Counts)

    ε = .8C.ll + rand()*.4C.ll # base step size ± 20%.
    loc_init(Y,S,X)
    loc_unconstrained(C,X)
    loc_momenta(C,X)
    ep = loc_potential(C,X) # previous energy; also computes the gradient.
    kp = loc_kinetic(X)

    # loc_println(X) # debug.
    # println("ep: $ep, kp: $kp") # debug.

    loc_pstep(.5ε,X)
    for i = 1:C.Ll
        loc_qstep(ε,Y,C,X)
        if i != C.Ll
            loc_potential(C,X)
            loc_pstep(ε,X)
        end
    end
    en = loc_potential(C,X)
    loc_pstep(.5ε,X)
    kn = loc_kinetic(X)

    # loc_println(X) # debug.
    # println("en: $en, kn: $kn") # debug.
    # println("diff: $(ep+kp-en-kn)") # debug.
    
    p = rand()
    l = exp(ep+kp-en-kn)
    # println((p < l) ? "$p < $l" : "$p ≥ $l") # debug.

    if p < l
        S.locs[:] = X.locs
        S.hts[:] = X.hts # order might have changed.
        U.laccepts += 1
    else
        loc_init(Y,S,X)        # leave the auxiliary parameters
        loc_unconstrained(C,X) # in a consistent state.
    end
    U.lprops += 1

end

function loc_pstep(s, X::Auxiliary)
    @. X.ptlocs -= s * X.∇tlocs
end

function loc_qstep(s, Y, C::Constants, X::Auxiliary)
    @. X.tlocs += s * X.ptlocs
    if !issorted(X.tlocs)
        sortperm!(X.ix, X.tlocs, initialized=true)
        X.tlocs[:] = X.tlocs[X.ix]
        X.ptlocs[:] = X.ptlocs[X.ix]
        X.hts[2:end] = @view(X.hts[2:end])[X.ix]
        X.lhts[2:end] = @view(X.lhts[2:end])[X.ix]
        X.ix[:] = 1:length(X.ix)
    end
    loc_constrained(Y,C,X)
end

#= Energy =====================================================================#

function loc_potential(C::Constants, X::Auxiliary)
    k,R,α,β = C.k, C.R, C.α, C.β
    hts, widths, counts = X.hts, X.widths, X.counts
    locs = @view(X.locs[2:end-1])

    e = sum(@. (widths+β)*hts - (counts+α)*X.lhts - log(widths)) -
        sum(@. log(locs*(R-locs)))
    X.∇tlocs .= @. begin
        (@view(hts[1:end-1]) - @view(hts[2:end]) -
        1/@view(widths[1:end-1]) + 1/@view(widths[2:end]))*(locs*(R-locs)/R) -
        (1-2locs/R)
    end
    return e
end

loc_kinetic(X::Auxiliary) = .5*sum(X.ptlocs.^2)

#= I/O ========================================================================#

function loc_println(X::Auxiliary)
    println("X: locs = $(string(X.locs)), hts = $(string(X.hts))\n",
            "   widths = $(string(X.widths)), counts = $(X.counts)\n",
            "   tlocs = $(string(X.tlocs)), ∇tlocs = $(string(X.∇tlocs))\n",
            "   ptlocs = $(string(X.ptlocs)), ix = $(X.ix)")
end