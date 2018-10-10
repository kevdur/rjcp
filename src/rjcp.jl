#===============================================================================
# Reversible Jump Change-point Analysis
# July 2018
#
# A relatively straightforward port of Green's Fortran implementation, but
# making use of Hamiltonian Monte Carlo for height updates. (One can also use
# HMC for location updates, but the discontinuities in the energy/likelihood
# function make it less appopriate.)
#
# The simplest way of running this program is to call `rjcp` with the output of
# `readvec`.
===============================================================================#

module RJCP

export green, readvec, rjcp

import Base: show, string

include("structs.jl") # 
include("k.jl")       # Must be included
include("loc.jl")     # in this order.
include("ht.jl")      #

#==== Main ====================================================================#

"""
    rjcp(filename; <keyword arguments>)

Performs reversible jump and Hamiltonian MCMC for the Poisson-process
change-point problem.

# Arguments
  - `data::Vector{Float64}`: an array of non-negative event times.
  - `R::Float64=40907.0`: an upper limit on the event times present in the data.
        The interval ``[0,R]`` will be used as an overall time window in the
        likelihood calculations.
  - `k::Union{Int,UnitRange{Int}}=0:20`: allowable values for the number of
        change points. If an integer is given, between-model jumps will not
        occur.
  - `loc::Symbol=:gibbs`, `ht::Symbol=:hmc`: the algorithms to be used for
        location and height updates respectively: `:gibbs` or `:hmc`.
  - `lambda::Float64=3.0`: the rate of the Poisson prior on the number of change
        points.
  - `alpha::Float64=1.0`, `beta::Float64=200.0`: the rates of the Beta prior on
        step heights.
  - `leps::Float64=.005`, `lleaps::Int=5`: leapfrog step size and count for
        location moves (when `loc` is set to `:hmc`).
  - `heps::Float64=.1`, `hleaps::Int=100`: leapfrog step size and count for
        height moves.
  - `samples::Int=40_000`, `burnin::Int=samples/10`: the number of samples to
        draw. Unless otherwise specified, a warm-up or 'burn-in' stage will be
        run prior to the sampling.
  - `frequency::Int=1`: the frequency at which samples are written to file,
        e.g., a frequency of 2 implies that every second sample is discarded.
"""
function rjcp(data::Vector{Float64};
              R::Float64=40907., k::Union{Int,UnitRange{Int}}=0:20,
              loc::Symbol=:gibbs, ht::Symbol=:hmc,
              lambda::Float64=3., alpha::Float64=1., beta::Float64=200.,
              leps::Float64=.005, lleaps::Int=5,
              heps::Float64=.1, hleaps::Int=100,
              samples::Int=40_000, burnin::Int=div(samples,10),
              frequency::Int=1)

    kmin,kmax = (k isa Int) ? (k,k) : extrema(k)
    Y = sort(data)
    C = Constants(kmax,R,lambda,alpha,beta,leps,lleaps,heps,hleaps)
    S = State(Y,kmin,kmax,R)

    X = Auxiliary(kmax)
    aux_init(Y,C,S,X)

    J = Jumps(k -> min(1,poissonratio(lambda,k)),
              k -> min(1,rpoissonratio(lambda,k-1)), kmin, kmax, .9)

    function rjmcmc(name::AbstractString, I::Int)
        (I <= 0) && return

        U = Counts(kmax)
        open(name*".txt", "w") do f
            write(f, header(kmax)*"\n")
            for i = 1:I
                if !isa(k,Int) rjmove(Y,C,S,X,J,U) end
                locmove(loc,Y,C,S,X,U)
                htmove(ht,Y,C,S,X,U)
                U.ks[S.k+1] += 1
                if (i % frequency == 0) write(f, string(S)*"\n") end
            end
        end
        println(U)
    end

    rjmcmc("burnin", burnin)
    rjmcmc("samples", samples)

end

function aux_init(Y,C,S,X)
    loc_init(Y,S,X)
    loc_unconstrained(S.k,C,X)
    ht_init(Y,S,X)
    ht_unconstrained(S.k,X)
end

#==== I/O =====================================================================#

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

"""
    green(name)

Converts output from Green's change-point program to the format of this program.

# Arguments
  - `name::AbstractString`: the shared of the two input files (`name.pos` and
        `name.ht`).
"""
function green(name::AbstractString)

    open(name*".pos") do pos
        open(name*".ht") do ht
            open(name*".txt", "w") do out

                values(l) = map(s -> parse(Float64,s), split(l))

                ps = [values(l) for l in readlines(pos)]
                hs = [values(l) for l in readlines(ht)]
                K = mapreduce(length,max,0,ps)
                S = State(readvec("coal.csv"),1,K,40907.)
                
                write(out, header(K)*"\n")
                for (i,p) in enumerate(ps)
                    k = length(p)
                    S.k = k
                    S.locs[1:k+2] .= vcat(0.,p,40907.)
                    S.hts[1:k+1] .= hs[i]
                    write(out, string(S)*"\n")
                end

            end
        end
    end

end

end