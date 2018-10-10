#=
=
= Filler.
=
=#

include("hcp.jl")

using StatPlots
using HCP: Constants, params_metrics,
           pstep, qstep, potential, kinetic, readvec

const Y = sort(readvec("coal.csv"))
const C = Constants(3,40907.,1.,200.,0.,0.)
const P = HCP.Parameters(3)
const n = 10_000
const xs = Vector{Float64}(n)
const es = Vector{Float64}(n)
const sh = [-1,-1]
const xd = [0.,0.]

const ls = :dash
const shape = :xcross
const ms = 5
const msw = 2

"""
Viewable parameters are s_i, i ∈ {1,2,3}, and h_j, j ∈ {0,1,2,3}.
"""
function start(sym=:s1)
    
    P.locs[1], P.locs[5] = 0., 40907.
    P.locs[2] = 14_500 + 1000randn()
    P.locs[3] = 24_000 + 3000randn()
    P.locs[4] = 35_500 + 1000randn()
    P.heights[1] = abs(.0085 + .001randn())
    P.heights[2] = abs(.0025 + .0005randn())
    P.heights[3] = abs(.003 + .001randn())
    P.heights[4] = abs(.001 + .0005randn())
    params_metrics(Y,P)
    potential(Y,C,P)

    println("s: $(P.locs)")
    println("h: $(P.heights)")

    if sym in [:s1,:s2,:s3]
        sh[1] = s = (sym == :s1) ? 1 :
                    (sym == :s2) ? 2 : 3
        xs[:] = collect(linspace(P.locs[s],P.locs[s+2],n))
        es[:] = eloc.(xs,s)
        xd = [P.locs[1+1], P.∇locs[s]]
        println("e: $(eloc(xd[1],s))")
    else
        sh[2] = h = (sym == :h0) ? 0 :
                    (sym == :h1) ? 1 :
                    (sym == :h2) ? 2 : 3
        xs[:] = linspace(.0001,.011,n)
        es[:] = eht.(xs,h)
        xd = [P.heights[h+1], P.∇heights[h+1]]
        println("e: $(eht(xd[1],h))")
    end
    show()

end

function step(;eps=1,N=100)

    ε = .8eps + rand()*.4eps # random interval ε ± 20%.

    if sh[1] > 0
        s = sh[1]
        
        P.plocs[s] = 1e-5randn()
        P.plocs[s] -= .5ε * P.∇locs[s] # pstep.
        println("p: $(P.plocs[s])")
        for i = 1:N
            P.locs[s+1] += ε * P.plocs[s] # qstep.
            params_metrics(Y,P)
            if i != N
                potential(Y,C,P)
                P.plocs[s] -= ε * P.∇locs[s] # pstep.
            end
        end
        potential(Y,C,P)

        xd = [P.locs[s+1], P.∇locs[s]]
        println("s: $(xd[1])")
    else
        h = sh[2] + 1

        P.pheights[h] = 10randn()
        P.pheights[h] -= .5ε * P.∇heights[h] # pstep.
        println("p: $(P.pheights[h])")
        for i = 1:N
            P.heights[h] += ε * P.pheights[h] # qstep.
            params_metrics(Y,P)
            if i != N
                potential(Y,C,P)
                P.pheights[h] -= ε * P.∇heights[h] # pstep.
            end
        end
        potential(Y,C,P)

        xd = [P.heights[h], P.∇heights[h]]
        println("h: $(xd[1])")
    end
    show()

end

function show()
    x,d = xd
    e = (sh[1] > 0) ? eloc(x,sh[1]) : eht(x,sh[2])
    m = searchsortedfirst(xs,x)
    i,j = clamp(m-n÷10,1,n), clamp(m+n÷10,1,n)
    plot(xs[i:j], es[i:j], leg=false, size=(1280,720))
    scatter!((x,e), leg=false, shape=shape, ms=ms, msw=msw)
    plot!(linspace(x-1000,x+1000,1000), z -> e+(z-x)d, leg=false, ls=ls)
end

function eloc(x,s)
    px = P.locs[s+1]
    P.locs[s+1] = (x < P.locs[s]) ? P.locs[s] :
                  (x > P.locs[s+2]) ? P.locs[s+2] : x
    params_metrics(Y,P)
    e = sum(@. begin
        (P.widths+C.β)*P.heights - (P.counts+C.α-1)*log(P.heights) -
        log(P.widths)
    end)
    P.locs[s+1] = px
    return e
end

function dloc(x,s)
    px = P.locs[s+1]
    P.locs[s+1] = (x < P.locs[s]) ? P.locs[s] :
                  (x > P.locs[s+2]) ? P.locs[s+2] : x
    params_metrics(Y,P)
    P.∇locs .= @. begin
        @view(P.heights[1:end-1]) - @view(P.heights[2:end]) -
        1/@view(P.widths[1:end-1]) + 1/@view(P.widths[2:end])
    end
    P.locs[s+1] = px
    return P.∇locs[s]
end

function eht(x,h)
    px = P.heights[h+1]
    P.heights[h+1] = x
    params_metrics(Y,P)
    e = sum(@. begin
        (P.widths+C.β)*P.heights - (P.counts+C.α-1)*log(P.heights) -
        log(P.widths)
    end)
    P.height[h+1] = px
    return e
end

function dht(x,h)
    px = P.heights[h+1]
    P.heights[h+1] = x
    params_metrics(Y,P)
    @. P.∇heights = (P.widths+C.β) - (P.counts+C.α-1)/P.heights
    P.heights[h+1] = px
    return P.∇heights[h+1]
end