
#==== Move ====================================================================#

function locmove(l::Symbol, Y, C::Constants, S::State, X::Auxiliary, U::Counts)
    if l == :hmc
        loc_hmc(Y,C,S,X,U)
    else
        loc_gibbs(Y,C,S,X,U)
    end
end

#==== Gibbs ===================================================================#

function loc_gibbs(Y, C::Constants, S::State, X::Auxiliary, U::Counts)

    k,α,β = S.k, C.α, C.β
    @views locs,hts =S.locs[1:k+2], S.hts[1:k+1]
    if (k == 0)
        return
    end

    @views ix = X.ix[2:k+1] # locations will be
    shuffle!(ix)            # updated in random order.
    for j in ix
        sp,s,sn = locs[j-1], locs[j], locs[j+1]
        hp,h = hts[j-1], hts[j]
        x = sp + rand()*(sn-sp)

        lratio = (x < s) ? llhdratio(Y,x,s,s,h,h,hp) :
                           llhdratio(Y,s,x,x,hp,hp,h)
        lratio += log((sn-x)/(sn-s) * (x-sp)/(s-sp))

        ratio = exp(max(-30,min(0,lratio)))
        if rand() < ratio
            locs[j] = x
            U.laccepts += 1
        end
        U.lprops += 1
    end
    loc_init(Y,S,X)          # leave the auxiliary
    loc_unconstrained(k,C,X) # variables in a 
    ix .= 2:k+1              # consistent state.

end

#==== HMC =====================================================================#

function loc_hmc(Y, C::Constants, S::State, X::Auxiliary, U::Counts)

    k = S.k
    ε = .8C.leps + rand()*.4C.leps # base step size ± 20%.
    loc_init(Y,S,X)
    loc_unconstrained(k,C,X)
    loc_momenta(k,X)
    ep = loc_potential(k,C,X) # previous energy; also computes the gradient.
    kp = loc_kinetic(k,X)

    loc_pstep(.5ε,k,X)
    for i = 1:C.lleaps
        loc_qstep(ε,k,Y,C,X,U)
        if i != C.lleaps
            loc_potential(k,C,X)
            loc_pstep(ε,k,X)
        end
    end
    en = loc_potential(k,C,X)
    loc_pstep(.5ε,k,X)
    kn = loc_kinetic(k,X)

    p = rand()
    l = exp(ep+kp-en-kn)

    if p < l
        @views S.locs[1:k+2] .= X.locs[1:k+2]
        @views S.hts[1:k+1] .= X.hts[1:k+1] # order might have changed.
        U.laccepts += 1
    else
        aux_init(Y,C,S,X) # leave the aux. variables in a consistent state.
    end
    U.lprops += 1

end

function loc_pstep(s, k, X::Auxiliary)
    @. @views X.ptlocs[1:k] -= s * X.∇tlocs[1:k]
end

function loc_qstep(s, k, Y, C::Constants, X::Auxiliary, U::Counts)
    @views tlocs,ptlocs = X.tlocs[1:k], X.ptlocs[1:k]
    @. tlocs += s*ptlocs
    if !issorted(tlocs)
        @views ix = X.ix[1:k]
        sortperm!(ix, tlocs, initialized=true)
        permute!(tlocs, ix)
        permute!(ptlocs, ix)
        @views permute!(X.hts[2:k+1], ix)
        @views permute!(X.lhts[2:k+1], ix)
        ix .= 1:k
        U.lsorts += 1
    end
    loc_constrained(k,Y,C,X)
end

function loc_potential(k, C::Constants, X::Auxiliary)
    R,α,β = C.R, C.α, C.β
    @views locs,hts,lhts = X.locs[2:k+1], X.hts[1:k+1], X.lhts[1:k+1]
    @views widths,counts = X.widths[1:k+1], X.counts[1:k+1]

    e = sum(@. (widths+β)*hts - (counts+α)*lhts - log(widths)) -
        sum(@. log(locs*(R-locs)))
    X.∇tlocs[1:k] .= @. begin
        (@view(hts[1:end-1]) - @view(hts[2:end]) -
        1/@view(widths[1:end-1]) + 1/@view(widths[2:end]))*(locs*(R-locs)/R) -
        (1-2locs/R)
    end
    return e
end

loc_kinetic(k, X::Auxiliary) = @views .5 * sum(X.ptlocs[1:k].^2)

#==== Auxiliary ===============================================================#

function loc_init(Y, S::State, X::Auxiliary) # locs, widths, counts.
    k = S.k
    @views X.locs[1:k+2] .= S.locs[1:k+2]
    loc_metrics(k,Y,X)
end

function loc_metrics(k, Y, X::Auxiliary) # widths, counts.
    sp = 0.
    y = 1 # index of the first uncounted datum.
    for i = 1:k+1
        s = X.locs[i+1]
        X.widths[i] = s-sp
        n = searchsortedfirst(@view(Y[y:end]), s) - 1
        X.counts[i] = n
        sp,y = s,y+n
    end
end

function loc_unconstrained(k, C::Constants, X::Auxiliary) # tlocs, texps, dlocs.
    @views locs = X.locs[2:k+1]
    @. X.tlocs[1:k] = log(locs/(C.R-locs))
end

function loc_constrained(k, Y, C::Constants, X::Auxiliary) # locs, metrics.
    @. @views X.locs[2:k+1] = C.R/(1+exp(-X.tlocs[1:k]))
    loc_metrics(k,Y,X)
end

function loc_momenta(k, X::Auxiliary) # ptlocs.
    @views randn!(X.ptlocs[1:k])
end

#==== I/O =====================================================================#

function loc_println(k, X::Auxiliary)
    println("X: locs = $(string(X.locs[1:k+2])),",
              " hts = $(string(X.hts[1:k+1]))\n",
            "   widths = $(string(X.widths[1:k+1])),",
              " counts = $(X.counts[1:k+1])\n",
            "   tlocs = $(string(X.tlocs[1:k])),",
              " ∇tlocs = $(string(X.∇tlocs[1:k]))\n",
            "   ptlocs = $(string(X.ptlocs[1:k])), ix = $(X.ix[1:k])")
end