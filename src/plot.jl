
#==== Plotting functions ======================================================#

using CSV
using DataFrames
using LaTeXStrings
using Query
using StatPlots

# function r(; name::AbstractString="samples", ks=0:Inf)

#     df = read(name*"txt")
    

# end

function ks(name::AbstractString="samples.txt")

    # We add the 'rows_for_type_detect' argument because there seems to be a bug
    # in the 'read' method: even with 'allowmissing=:all', a column whose first
    # missing value only appears after the type-detect rows will cause an error.

    df = CSV.read(name, rows_for_type_detect=typemax(Int))
    N = size(df,1)
    df = by(df, :k, d -> DataFrame(p=size(d,1)/N))
    @df df bar(:k, :p, xlabel="k", ylabel="P(k)", label="", size=(600,400))

end

function locs(name::AbstractString="samples.txt"; k::Int=1, ss=1:k)

    df = CSV.read(name, rows_for_type_detect=typemax(Int))
    df = @from r in df begin
         @where r.k == k
         @select r
         @collect DataFrame
    end
    ss = [Symbol("s$i") for i in ss]
    ps = vcat(
        [density(df[s], title=s, label="") for s in ss],
        [scatter(df[s], label="", ms=1) for s in ss]
    )
    m = length(ss)
    plot(ps..., layout=(2,m), size=(600m, 600))

end

function hts(name::AbstractString="samples.txt"; k::Int=1, hs=0:k)

    df = CSV.read(name, rows_for_type_detect=typemax(Int))
    df = @from r in df begin
         @where r.k == k
         @select r
         @collect DataFrame
    end
    hs = [Symbol("h$i") for i in hs]
    ps = vcat(
        [density(df[h], title=h, label="") for h in hs],
        [scatter(df[h], label="", ms=1) for h in hs]
    )
    m = length(hs)
    plot(ps..., layout=(2,m), size=(600m, 600))

end

# using Plots

# import Base.linspace

# #= I/O functions ==============================================================#

# function load(name::AbstractString)

#     ks = Vector{Int}()
#     open(name*".k") do f
#         ks = eval(parse(readline(f)))
#     end

#     K = maximum(find(ks)) - 1 # maximum sampled k.
#     S = [Matrix{Float64}(k+2,ks[k+1]) for k = 0:K] # step locations.
#     H = [Matrix{Float64}(k+1,ks[k+1]) for k = 0:K] # step heights.
#     js = zeros(Int, K+1) # track the samples per k.
#     open(name*".txt") do f
#         for (l,ln) in enumerate(eachline(f))
#             s,h = map(t -> eval(parse(t)), split(ln,';'))
#             k = length(s) - 2
#             j = (js[k+1] += 1)
#             S[k+1][:,j] = s
#             H[k+1][:,j] = h
#         end
#     end
#     return js,S,H

# end

# #= Statistical functions ======================================================#

# # Returns a range of linearly spaced elements over the data domain.
# linspace(n::Int, S) = linspace(0, S[end][end,end]-1, n)

# # Returns the corresponding step heights for the given (sorted) datums.
# function heights(xs,ss,hs)
#     is = similar(xs)
#     c = 1 # cursor.
#     for (i,s) in enumerate(@view(ss[2:end]))
#         j = c-1 + searchsortedfirst(@view(xs[c:end]), s)
#         is[c:j-1] = fill(hs[i],j-c)
#         c = j
#     end
#     return is
# end

# # Returns the sum of the step heights for given data, for a specific k.
# function hsum(xs,k,S,H)
#     ss,hs = S[k+1], H[k+1]
#     return (size(ss,2) == 0) ?
#         zeros(xs) : sum(i -> heights(xs,ss[:,i],hs[:,i]), indices(ss,2))
# end

# # Returns the posterior rate function at the given datums, for a specific k.
# function r(xs,k,S,H)
#     kn = size(S[k+1],2)
#     return (kn == 0) ? zeros(xs) : hsum(xs,k,S,H) / size(S[k+1],2)
# end

# # Returns the posterior rate function at the given datums.
# function r(xs,S,H)
#     K = length(S) - 1
#     ks = [size(S[k+1],2) for k = 0:K]
#     sum(k -> hsum(xs,k,S,H), 0:K) / sum(ks)
# end

# #= Plotting functions =========================================================#

# plot(xs, S, H) = Plots.plot(xs, r(xs,S,H), leg=false)
# plot(xs, k::Int, S, H) = Plots.plot(xs, r(xs,k,S,H), lab="k = $k")
# plot(xs, ks, S, H) = Plots.plot(xs, [r(xs,k,S,H) for k in ks],
#                                 lab=["k = $k" for k in ks])