
#==== Moves ===================================================================#

function rjmove(Y, C::Constants, S::State, X::Auxiliary, J::Jumps, U::Counts)
    k,bprob,dprob = S.k, J.bprob, J.dprob

    p = rand()
    b = bprob(k)
    if p < b
        bmove(Y,C,S,X,J,U)
    elseif p < b + dprob(k)
        dmove(Y,C,S,X,J,U)
    end
end

function bmove(Y, C::Constants, S::State, X::Auxiliary, J::Jumps, U::Counts)
    k,R,λ,α,β,lλ,αlβ = S.k, C.R, C.λ, C.α, C.β, C.lλ, C.αlβ

    s = rand()*R
    @views j = searchsortedlast(S.locs[1:k+2],s)
    s1,s2 = S.locs[j], S.locs[j+1]
    w,w1,w2 = s2-s1, s-s1, s2-s
    h,u = S.hts[j], rand()
    h1 = h * (u/(1-u))^(w2/w)
    h2 = h1 * ((1-u)/u)

    lratio = llhdratio(Y,s1,s2,s,h1,h2,h) +
             lpriorratio(k,h1,h2,h,w1,w2,w,C) +
             lproposalratio(J.bprob,J.dprob,k,C.R) +
             ljacobian(h1,h2,h)
    ratio = exp(clamp(lratio,-30,0))

    if rand() < ratio
        S.locs[k+3] = S.locs[k+2]
        for i = k+2:-1:j+2
            S.locs[i] = S.locs[i-1] # this can't be done
            S.hts[i] = S.hts[i-1]   # using .=, as in 'dmove'.
        end
        S.hts[j] = h1
        S.hts[j+1] = h2
        S.locs[j+1] = s
        S.k += 1
        U.baccepts += 1
        aux_init(Y,C,S,X) # leave the aux. variables in a consistent state.
    end
    U.bprops += 1
end

function dmove(Y, C::Constants, S::State, X::Auxiliary, J::Jumps, U::Counts)
    k,R,λ,α,β,lλ,αlβ = S.k, C.R, C.λ, C.α, C.β, C.lλ, C.αlβ

    j = rand(2:k+1)
    s1,s,s2 = S.locs[j-1], S.locs[j], S.locs[j+1]
    w,w1,w2 = s2-s1, s-s1, s2-s
    h1,h2 = S.hts[j-1], S.hts[j]
    h = h1^(w1/w) * h2^(w2/w)

    lratio = llhdratio(Y,s1,s2,s,h1,h2,h) +
             lpriorratio(k-1,h1,h2,h,w1,w2,w,C) +
             lproposalratio(J.bprob,J.dprob,k,C.R) +
             ljacobian(h1,h2,h)
    ratio = exp(clamp(-lratio,-30,0)) # note the negation.

    if rand() < ratio
        @views S.locs[j:k+1] .= S.locs[j+1:k+2]
        @views S.hts[j:k] .= S.hts[j+1:k+1]
        S.locs[k+2] = NaN
        S.hts[k+1] = NaN
        S.hts[j-1] = h
        S.k -= 1
        U.daccepts += 1
        aux_init(Y,C,S,X)
    end
    U.dprops += 1
end

#==== Probability ratios ======================================================#

"""
Returns the ratio p(k+1)/p(k) of two Poisson(λ) probabilities.
"""
poissonratio(λ,k) = λ/(k+1)
rpoissonratio(λ,k) = (k+1)/λ
lpoissonratio(lλ,k) = lλ - log(k+1)

"""
Returns the ratio p(w1,w2,...|k+1)/p(w,...|k) of the probabilities of two
even-numbered order statistics.
"""
lstatsratio(R,k,w1,w2,w) = log((2k+2)*(2k+3)/R^2 * w1*w2/w)

"""
Returns the natural logarithm of the ratio p(h1,h2,...|k+1)/p(h,...|k) of two
Beta(α,β) probabilities.
"""
lbetaratio(α,β,αlβ,h1,h2,h) = αlβ + (α-1)*log(h1*h2/h) - (h1+h2-h)β

"""
Returns the natural log-likelihood ratio required for a birth move.

# Arguments
  - `s1`, `s2`, `s`: two existing, ordered step locations and a new one that
        lies between them.
  - `h1`, `h2`, `h`: the heights of the new steps adjacent to `s`, as well as
        the height of the step that they are replacing.
"""
function llhdratio(Y,s1,s2,s,h1,h2,h)
    n1,n2,n3 = [searchsortedfirst(Y,s) for s in (s1,s,s2)]
    return (n2-n1)*log(h1/h) - (h1-h)*(s-s1) +
           (n3-n2)*log(h2/h) - (h2-h)*(s2-s)
end

"""
Returns the natural logarithm of the prior ratio for a birth move.

# Arguments
  - `h1`, `h2`, `h`: the heights of two new adjacent steps, along with the
        height `h` of the step that they are replacing.
  - `w1`, `w2`, `w`: the widths of the steps of heights `h1`, `h2`, and `h`,
        respectively.
"""
function lpriorratio(k,h1,h2,h,w1,w2,w,C::Constants)
    return lpoissonratio(C.lλ,k) +
           lstatsratio(C.R,k,w1,w2,w) +
           lbetaratio(C.α,C.β,C.αlβ,h1,h2,h)
end

"""
Returns the natural logarithm of the proposal ratio for a birth move.
"""
lproposalratio(b,d,k,R) = log(d(k+1) * R/(b(k)*(k+1)))

"""
Returns the natural logarithm of the Jacobian factor for a birth move.
"""
ljacobian(h1,h2,h) = log((h1+h2)^2 / h)

#==== I/O =====================================================================#

function k_println(S::State)
    println("S: k = $S.k,",
              " locs = $(string(S.locs[1:k+2])),",
              " hts = $(string(S.hts[1:k+1]))\n")
end

function k_println(k, X::Auxiliary)
    println("X: locs = $(string(X.locs[1:k+2])),",
              " hts = $(string(X.hts[1:k+1]))\n",
            "   widths = $(string(X.widths[1:k+1])),",
              " counts = $(X.counts[1:k+1])\n",
            "   tlocs = $(string(X.tlocs[1:k])),",
              " lhts = $(string(X.lhts[1:k+1]))\n",
            "   ix = $(X.ix[1:k])\n")
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