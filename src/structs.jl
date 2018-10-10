
#==== Command line arguments ==================================================#

struct Constants
    K::Int
    R::Float64
    λ::Float64
    α::Float64
    β::Float64
    leps::Float64
    lleaps::Int
    heps::Float64
    hleaps::Int
    lλ::Float64  # used in prior
    αlβ::Float64 # calculations.
    Constants(K,R,λ,α,β,leps,heps,lleaps,hleaps) = new(
        K,R,λ,α,β,leps,heps,lleaps,hleaps,
        log(λ),
        α*log(β)-lgamma(α)
    )
end

#==== State ===================================================================#

mutable struct State
    k::Int
    locs::Vector{Float64}
    hts::Vector{Float64}
end

"""
Initialises a new state with evenly spaced step locations and mean heights.
"""
function State(Y::Vector{Float64}, k::Int, K::Int, R::Float64)
    S = State(k, fill(NaN,K+2), fill(NaN,K+1))
    S.locs[1:k+2] .= linspace(0.,R,k+2)
    s,y = 0.,1 # step location and data index.
    for i = 1:k+1
        sn = S.locs[i+1] # next step location.
        @views n = searchsortedfirst(Y[y:end], sn) - 1 # step data count.
        S.hts[i] = n/(sn-s)
        s,y = sn, y+n
    end
    return S
end

function header(K)
    io = IOBuffer()
    write(io, "k,")
    for k = 0:K+1
        write(io, "s$k,")
    end
    for k = 0:K
        write(io, "h$k")
        if (k < K) write(io, ",") end
    end
    return String(take!(io))
end

function string(S::State)
    k = S.k
    io = IOBuffer()
    write(io, string(k,","))
    for s in S.locs
        if !isnan(s) write(io, string(s)) end
        write(io, ",")
    end
    for (i,h) in enumerate(S.hts)
        if !isnan(h) write(io, string(h)) end
        if (i != length(S.hts)) write(io, ",") end
    end
    return String(take!(io))
end

#==== Extended state ==========================================================#

struct Auxiliary
    locs::Vector{Float64}
    hts::Vector{Float64}
    widths::Vector{Float64}
    counts::Vector{Int}
    tlocs::Vector{Float64} # t = logit(s/R).
    lhts::Vector{Float64}  # j = ln(h).
    ∇tlocs::Vector{Float64}
    ∇lhts::Vector{Float64}
    ptlocs::Vector{Float64}
    plhts::Vector{Float64}
    ix::Vector{Int} # index vector used for sorting.
    Auxiliary(K) = new(
        Vector{Float64}(K+2),
        Vector{Float64}(K+1),
        Vector{Float64}(K+1),
        Vector{Int}(K+1),
        Vector{Float64}(K),
        Vector{Float64}(K+1),
        Vector{Float64}(K),
        Vector{Float64}(K+1),
        Vector{Float64}(K),
        Vector{Float64}(K+1),
        collect(1:K+1)
    )
end

#==== Birth/death probabilities ===============================================#

struct Jumps
    bprob::Function
    dprob::Function
end

"""
Initialises a new pair of jump proposals such that the sum of their returned
probabilities never exceeds the given value.
"""
function Jumps(b::Function, d::Function, kmin::Int, kmax::Int, s::Float64=.9)
    c = max(b(kmin),d(kmax))
    if kmin < kmax-1
        c = max(c, maximum(b(k)+d(k) for k = kmin+1:kmax-1))
    end
    c = s/c
    Jumps(k -> (k < kmax) ? c*b(k) : 0,
          k -> (k > kmin) ? c*d(k) : 0)
end

#==== Move statistics =========================================================#

mutable struct Counts
    bprops::Int
    dprops::Int
    lprops::Int
    hprops::Int
    baccepts::Int
    daccepts::Int
    laccepts::Int
    haccepts::Int
    lsorts::Int
    ks::Vector{Int}
    Counts(K) = new(0,0,0,0,0,0,0,0,0,zeros(Int,K+1))
end

function show(io::IO, C::Counts)

    K = length(C.ks) - 1
    println(io, "No. of change points:")
    for k = 0:K
        m = C.ks[k+1]
        (m > 0) && println(io, "  $k: $m")
    end
    println(io, "Mean: $(dot(C.ks,0:K) / sum(C.ks))")

    println(io, "Acceptance rates:")
    @printf(io, "  birth: %d/%d (%.2f%%)\n",
        C.baccepts, C.bprops, (C.bprops == 0) ? 0 : 100*C.baccepts/C.bprops)
    @printf(io, "  death: %d/%d (%.2f%%)\n",
        C.daccepts, C.dprops, (C.dprops == 0) ? 0 : 100*C.daccepts/C.dprops)
    @printf(io, "  location: %d/%d (%.2f%%)\n",
        C.laccepts, C.lprops, (C.lprops == 0) ? 0 : 100*C.laccepts/C.lprops)
    @printf(io, "  height: %d/%d (%.2f%%)",
        C.haccepts, C.hprops, (C.hprops == 0) ? 0 : 100*C.haccepts/C.hprops)
    if C.lsorts > 0
        print(io, "\nLocation sorts: $(C.lsorts)")
    end
    
end