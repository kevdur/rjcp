#=
= Simulation state structures and methods.
=#

#= Misc. ======================================================================#

mutable struct Counts
    lprops::Int
    hprops::Int
    laccepts::Int
    haccepts::Int
    Counts() = new(0,0,0,0)
end

struct Constants
    k::Int
    R::Float64
    α::Float64
    β::Float64
    ll::Float64 # leapfrog step sizes for
    lh::Float64 # location and height moves.
    Ll::Int     # leapfrog step counts for
    Lh::Int     # location and height moves.
end

#= State ======================================================================#

struct State
    locs::Vector{Float64}
    hts::Vector{Float64}
end

"""
Initialises a new state with evenly spaced step locations and average rates.
"""
function State(Y::Vector{Float64}, k::Int, R::Float64)
    S = State(collect(linspace(0.,R,k+2)), Vector{Float64}(k+1))
    s,y = 0., 1 # step location and data index.
    for i = 1:k+1
        sn = S.locs[i+1] # next step location.
        n = searchsortedfirst(@view(Y[y:end]), sn) - 1 # step data count.
        S.hts[i] = n/(sn-s)
        s,y = sn, y+n
    end
    return S
end

function string(S::State)
    io = IOBuffer()
    for s in S.locs
        write(io, string(s,", "))
    end
    for h in @view(S.hts[1:end-1])
        write(io, string(h,", "))
    end
    write(io, string(S.hts[end]))
    return String(take!(io))
end

#= Extended state =============================================================#

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
    Auxiliary(k) = new(
        Vector{Float64}(k+2),
        Vector{Float64}(k+1),
        Vector{Float64}(k+1),
        Vector{Int}(k+1),
        Vector{Float64}(k),
        Vector{Float64}(k+1),
        Vector{Float64}(k),
        Vector{Float64}(k+1),
        Vector{Float64}(k),
        Vector{Float64}(k+1),
        collect(1:k)
    )
end