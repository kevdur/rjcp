#=
= Energy and derivative plots.
=#

include("hcp.jl")

using HCP: Constants, Auxiliary, readvec
using LaTeXStrings
using StatPlots

function loc(data::Vector{Float64};
    s::Int=1, R::Float64=40907.,
    alpha::Float64=1., beta::Float64=200.)

    Y = sort(data)
    C = Constants(3,40907.,1.,200.,0.,0.,0.,0.)
    X = Auxiliary(3)
    X.locs[1], X.locs[5] = 0., 40907.
    X.locs[2] = 14_500 + 1000randn()
    X.locs[3] = 24_000 + 3000randn()
    X.locs[4] = 35_500 + 1000randn()
    X.hts[1] = abs(.0085 + .001randn())
    X.hts[2] = abs(.0025 + .0005randn())
    X.hts[3] = abs(.003 + .001randn())
    X.hts[4] = abs(.001 + .0005randn())
    HCP.loc_metrics(Y,X)
    HCP.loc_unconstrained(C,X)
    HCP.loc_potential(C,X)
    HCP.loc_momenta(C,X)
    HCP.ht_unconstrained(X)
    HCP.ht_momenta(C,X)
    HCP.loc_println(X)

    function e(x)
        X.locs[s+1] = x
        HCP.loc_metrics(Y,X)
        HCP.loc_unconstrained(C,X)
        return HCP.loc_potential(C,X)
    end

    function d(x)
        X.locs[s+1] = x
        HCP.loc_metrics(Y,X)
        HCP.loc_unconstrained(C,X)
        HCP.loc_potential(C,X)
        return X.∇tlocs[s]
    end

    r = X.locs[s+2]-X.locs[s]
    xp,xn = X.locs[s]+r/50, X.locs[s+2]-r/50
    ts = linspace(log(xp/(R-xp)), log(xn/(R-xn)), 1000)
    xs = map(t -> R/(1+exp(-t)), ts)
    es,ds = e.(xs), d.(xs)
    
    # Approximate the integral of the derivative to compare it to the
    # discontinuous function.

    fs = similar(es)
    fs[1] = 0
    for i = 2:length(ts)
        fs[i] = fs[i-1] + (ts[i]-ts[i-1])ds[i-1]
    end
    j = div(length(ts),2)
    fs += es[j]-fs[j]

    # Compute the derivatives numerically to verify the analyitic solution.
    # Sampled outliers arising from discontinuities in the energy function are
    # removed (hopefully).

    tmids = @. @view(ts[1:end-1]) + .5(@view(ts[2:end]) - @view(ts[1:end-1]))
    dmids = @. begin
        (@view(es[2:end]) - @view(es[1:end-1])) /
        (@view(ts[2:end]) - @view(ts[1:end-1]))
    end

    L,U = extrema(ds)
    for (i,d) in enumerate(dmids)
        (d < L || d > U) && (dmids[i] = NaN)
    end

    # Plot the numerical and analytical derivatives.

    p = plot(ts,es,lab="energy")
    plot!(ts,fs,lab="integral",ls=:dash)
    q = plot(ts,ds,lw=8,lab="analytic")
    plot!(q,tmids,dmids,lw=4,lab="numeric")
    plot(p,q,layout=(2,1),size=(1200,900))

    # p = plot(ts,es,lc=:black,lab=false,ylabel=L"H")
    # plot!(ts,fs,lc=:black,lab=false,ls=:dash)
    # q = plot(ts,ds,lc=:black,lab=false,xlabel=L"\hat{s}_2",ylabel=L"\frac{\partial H}{\partial \hat{s}_2}")
    # pq = plot(p,q,layout=(2,1),size=(480,360))
    # savefig(pq, "loc.pdf")
    # savefig(pq, "loc.png")
    # savefig(pq, "loc.tex")

end

function ht(data::Vector{Float64};
    h::Int=0, R::Float64=40907.,
    alpha::Float64=1., beta::Float64=200.)

    Y = sort(data)
    C = Constants(3,40907.,1.,200.,0.,0.,0.,0.)
    X = Auxiliary(3)
    X.locs[1], X.locs[5] = 0., 40907.
    X.locs[2] = 14_500 + 1000randn()
    X.locs[3] = 24_000 + 3000randn()
    X.locs[4] = 35_500 + 1000randn()
    X.hts[1] = abs(.0085 + .001randn())
    X.hts[2] = abs(.0025 + .0005randn())
    X.hts[3] = abs(.003 + .001randn())
    X.hts[4] = abs(.001 + .0005randn())
    HCP.loc_metrics(Y,X)
    HCP.loc_unconstrained(C,X)
    HCP.loc_potential(C,X)
    HCP.loc_momenta(C,X)
    HCP.ht_unconstrained(X)
    HCP.ht_momenta(C,X)
    HCP.ht_println(X)

    function e(x)
        X.hts[h+1] = x
        HCP.ht_unconstrained(X)
        return HCP.ht_potential(C,X)
    end

    function d(x)
        X.hts[h+1] = x
        HCP.ht_unconstrained(X)
        HCP.ht_potential(C,X)
        return X.∇lhts[h+1]
    end

    js = linspace(log(.0001),log(.025),1000)
    xs = exp.(js)
    es,ds = e.(xs), d.(xs)
    jmids = @. @view(js[1:end-1]) + .5(@view(js[2:end]) - @view(js[1:end-1]))
    dmids = @. begin
        (@view(es[2:end]) - @view(es[1:end-1])) /
        (@view(js[2:end]) - @view(js[1:end-1]))
    end

    p = plot(js,es,lab="energy")
    q = plot(js,ds,lw=8,lab="analytic")
    plot!(q,jmids,dmids,lw=4,lab="numeric")
    plot(p,q,layout=(2,1),size=(1200,900))

end