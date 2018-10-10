
#==== Move ====================================================================#

function htmove(h::Symbol, Y, C::Constants, S::State, X::Auxiliary, U::Counts)
    if h == :hmc
        ht_hmc(Y,C,S,X,U)
    else
        ht_gibbs(Y,C,S,X,U)
    end
end

#==== Gibbs ===================================================================#

function ht_gibbs(Y, C::Constants, S::State, X::Auxiliary, U::Counts)

    k,α,β = S.k, C.α, C.β
    @views locs,hts = S.locs[1:k+2], S.hts[1:k+1]

    @views ix = X.ix[1:k+1] # heights will be
    shuffle!(ix)            # updated in random order.
    for j in ix
        s,sn = locs[j], locs[j+1]
        h = hts[j]
        x = h*exp(rand()-.5)

        lratio = llhdratio(Y,s,sn,sn,x,x,h) + α*log(x/h) - (x-h)β

        ratio = exp(max(-30,min(0,lratio)))
        if rand() < ratio
            hts[j] = x
            U.haccepts += 1
        end
        U.hprops += 1
    end
    ht_init(Y,S,X)        # leave the auxiliary
    ht_unconstrained(k,X) # variables in a 
    ix .= 1:k+1           # consistent state.

end

#==== HMC =====================================================================#

function ht_hmc(Y, C::Constants, S::State, X::Auxiliary, U::Counts)

    k = S.k
    ε = .8C.heps + rand()*.4C.heps # base step size ± 20%.
    ht_init(Y,S,X)
    ht_unconstrained(k,X)
    ht_momenta(k,X)
    ep = ht_potential(k,C,X) # previous energy; also computes the gradient.
    kp = ht_kinetic(k,X)

    ht_pstep(.5ε,k,X)
    for i = 1:C.hleaps
        ht_qstep(ε,k,X)
        if i != C.hleaps
            ht_potential(k,C,X)
            ht_pstep(ε,k,X)
        end
    end
    en = ht_potential(k,C,X)
    ht_pstep(.5ε,k,X)
    kn = ht_kinetic(k,X)

    p = rand()
    l = exp(ep+kp-en-kn)

    if p < l
        @views S.hts[1:k+1] .= X.hts[1:k+1]
        U.haccepts += 1
    else
        ht_init(Y,S,X)        # leave the auxiliary variables
        ht_unconstrained(k,X) # in a consistent state.
    end
    U.hprops += 1

end

function ht_pstep(s, k, X::Auxiliary)
    @. @views X.plhts[1:k+1] -= s * X.∇lhts[1:k+1]
end

function ht_qstep(s, k, X::Auxiliary)
    @. @views X.lhts[1:k+1] += s * X.plhts[1:k+1]
    ht_constrained(k,X)
end

function ht_potential(k, C::Constants, X::Auxiliary)
    R,α,β = C.R, C.α, C.β
    @views locs,hts,lhts = X.locs[2:k+1], X.hts[1:k+1], X.lhts[1:k+1]
    @views widths,counts = X.widths[1:k+1], X.counts[1:k+1]

    e = sum(@. (widths+β)*hts - (counts+α)*lhts - log(widths)) -
        sum(@. log(locs*(R-locs)))
    @. X.∇lhts[1:k+1] = (widths+β)*hts - (counts+α)
    return e
end

ht_kinetic(k, X::Auxiliary) = @views .5 * sum(X.plhts[1:k+1].^2)

#==== Auxiliary ===============================================================#

function ht_init(Y, S::State, X::Auxiliary) # hts.
    k = S.k
    @views X.hts[1:k+1] .= S.hts[1:k+1]
end

function ht_unconstrained(k, X::Auxiliary) # lhts.
    @. @views X.lhts[1:k+1] = log(X.hts[1:k+1])
end

function ht_constrained(k, X::Auxiliary) # hts.
    @. @views X.hts[1:k+1] = exp(X.lhts[1:k+1])
end

function ht_momenta(k, X::Auxiliary) # plhts
    @views randn!(X.plhts[1:k+1])
end

#===== I/O =====================================================================#

function ht_println(k, X::Auxiliary)
    println("X: locs = $(string(X.locs[1:k+2])),",
              " hts = $(string(X.hts[1:k+1]))\n",
            "   lhts = $(string(X.lhts[1:k+1])),",
              " ∇lhts = $(string(X.∇lhts[1:k+1])),",
              " plhts = $(string(X.plhts[1:k+1]))")
end