#===============================================================================
# Hamiltonian Monte Carlo for the Poisson change-point problem with a fixed
# number of change points, k, and a step-like rate function.
#
# In this implementation, step locations and heights are transformed to
# unconstrained spaces, but the locations are not done so in an ordered
# manner -- their ordering must be checked after each leapfrog step.
#
# Because the locations and heights (and their derivatives) operate on different
# scales, the HMC move is split into two parts.
#
# Usage:
# To be honest, this implementation didn't perform as well as I initially
# thought it would -- Green's Gibbs-based implementation converges more quickly
# as far as I can tell. I believe the reason for this (which I only realised
# quite far along) is that the energy/likelihood function is riddled with
# discontinuities -- one for each datum -- and is almost step like in
# appearance. This implies that the gradient -- which is what provides HMC with
# its ability to avoid random walks -- doesn't provide information about the
# overall behaviour of the function, and is of relatively little use. It is for
# this reason that I've set the default number of leapfrog steps to one for the
# location move. The actual operational issues are a bit more subtle: the
# gradient usually informs the momentum variables, which compensate for the
# expected change in energy; in this case, these two changes disagree with each
# other across every discontinuity, and thus meaningful changes to the momentum
# often harm the performance of the simulation.
#
# In any case, the algorithm is still correct (as far as I can tell), and does
# work. A nice way to visualise its output is as follows (e.g., when k = 2):
#
#     using CSV,DataFrames,StatPlots
#     df = CSV.read("samples.txt", header=["s0","s1","s2","s3","h0","h1","h2"])
#     cs = (:s1,:s2) # or c = (:h0,:h1,:h2).
#     ps = vcat([density(hdf[c]) for c in cs], [scatter(hdf[c]) for c in cs])
#     plot(ps...)
===============================================================================#

module HCP

export hcp, readvec

import Base: show, string

include("structs.jl") # must appear
include("ht.jl")      # in this
include("loc.jl")     # order.

#= Main =======================================================================#

"""
    hcp(filename; <keyword arguments>)

Performs HMC for the Poisson-process change-point problem.

# Arguments
  - `data::Vector{Float64}`: an array of non-negative event times.
  - `k::Int=3`: the number of change points.
  - `R::Float64=40907.0`: an upper limit on the event times present in the data.
  - `alpha::Float64=1.0`, `beta::Float64=200.0`: the rates of the Beta prior on
        step heights.
  - `ll::Float64=.006`, `lh::Int=.1`: leapfrog step sizes for location and
        height moves respectively.
  - `Ll::Float64=1`, `Lh::Int=100`: leapfrog step counts for location and
        height moves respectively.
  - `samples::Int=10_000`, `burnin::Int=samples/10`: the number of samples to
        draw. Unless otherwise specified, a warm-up or 'burn-in' stage will be
        run prior to the sampling.
  - `frequency::Int=1`: the frequency at which samples are written to file,
        e.g., a frequency of 2 implies that every second sample is discarded.
"""
function hcp(data::Vector{Float64};
             k::Int=1, R::Float64=40907.,
             alpha::Float64=1., beta::Float64=200.,
             ll::Float64=.006, lh::Float64=.1,
             Ll::Int=1, Lh::Int=100,
             samples::Int=10_000, burnin::Int=div(samples,10),
             frequency::Int=1)

    Y = sort(data)
    C = Constants(k,R,alpha,beta,ll,lh,Ll,Lh)
    S = State(Y,k,R)
    
    X = Auxiliary(k)
    loc_init(Y,S,X)
    loc_unconstrained(C,X)
    ht_init(Y,S,X)
    ht_unconstrained(X)

    function hmc(name::AbstractString, I::Int)
        (I <= 0) && return

        U = Counts()
        open(name*".txt", "w") do f
            for i = 1:I
                loc_move(Y,C,S,X,U)
                ht_move(Y,C,S,X,U)
                (i % frequency == 0) && write(f, string(S)*"\n")
            end
        end
        println(name*":")
        println("  locs: $(U.laccepts)/$(U.lprops) " *
            @sprintf("(%.2f%%)", 100*U.laccepts/U.lprops))
        println("  hts:  $(U.haccepts)/$(U.hprops) " *
            @sprintf("(%.2f%%)", 100*U.haccepts/U.hprops))
    end

    hmc("burnin", burnin)
    hmc("samples", samples)

end 

#= I/O ========================================================================#

"""
    readvec(path)

Returns a `Vector{Float64}` of the numbers contained in the given file.

# Arguments
  - `path::Union{AbstractString,IO}`: the name or `IO` instance of the file
        containing the change-point data. These should be non-negative numbers
        formatted in a single row, column, or a matrix with no missing elements.
        Columns should be comma-separated.
"""
function readvec(path::Union{AbstractString,IO})
    vec(readdlm(path, ',', Float64))
end

end # module.