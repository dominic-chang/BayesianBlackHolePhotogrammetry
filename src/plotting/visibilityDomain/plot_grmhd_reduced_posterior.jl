using Revise
using Serialization, MCMCChains
using Comrade
using Pyehtim
using Plots
using Pigeons
using CairoMakie
using PairPlots
using ColorSchemes
using FileIO
using Krang
include(joinpath(dirname(@__DIR__) , "..", "models", "JuKeBOX.jl"))

red_cb = colorant"rgba(84%, 11%, 38%, 1.0)";
blue_cb = colorant"rgba(12%, 53%, 89%, 1.0)";
orange_cb = colorant"rgba(100%, 75%, 3%, 1.0)";
blue_cb_t = colorant"rgba(12%, 53%, 89%, 0.5)";
orange_cb_t = colorant"rgba(100%, 75%, 3%, 0.5)";
blue_cb_t1 = colorant"rgba(12%, 53%, 89%, 0.25)";
orange_cb_t1 = colorant"rgba(100%, 75%, 3%, 0.25)";
white_cb_t = colorant"rgba(100%, 100%, 100%, 0.5)";
green_cb = colorant"rgba(0%, 30%, 25%, 1.0)"

model_name = "sa+0.94_r160_GRMHD"
in_base = abspath(dirname(@__DIR__), "..", "visibilityDomain","results","all","2024-02-09-22-51-36-jUAZ7LdC")
include(joinpath(in_base, "model_params.jl"))
outpath = abspath(dirname(@__DIR__), ".." , "..","..","..","runs","visibilityDomain","data","$(model)")
#pt = PT(in_base)
#s_array = sample_array(pt)
rnd = 14
reduced_recorders = Serialization.deserialize(joinpath(in_base, "round=$rnd", "checkpoint", "reduced_recorders.jls"))
s_array = hcat([reduced_recorders.traces[n_tempering_levels=>i] for i in 1:2^rnd]...)'[:,1:end]

Plots.scatter(s_array[:, 1, :])
best_fit = begin 
    param_file = open(abspath(dirname(@__DIR__), "..", "..","runs","image_domain","sa+0.94_r160_nall_tavg.fits", "JBOX", "best_nxcorr.txt"))
    [readline(param_file) for _ in 1:3]
    eval(Meta.parse(readline(param_file)[12:end]))
end

nxcorr_vals = begin
    b_f_keys = keys(best_fit)[begin:end-2]
    NamedTuple{replace(b_f_keys, :σ=>:spec)}(map(x-> x ∈ (:θo, :θs, :χ, :ι, :η) ? best_fit[x]*180/π : best_fit[x], b_f_keys))
end
true_vals = (m_d=3.83, spin=-0.94, θo=17 / 180, pa=360 - 72)
prior_keys = collect(keys(prior))
samples = Chains(s_array, collect(prior_keys))
tsamples = Chains(reshape(hcat(collect.(map(x -> values(transform(cpost, x)), [samples.value[i, :] for i in 1:size(samples.value)[1]]))...)', size(samples.value)), collect(prior_keys))
MCMCChains.hpd(tsamples)

samples_to_plot= begin 
    temp = MCMCChains.to_matrix(tsamples[prior_keys])

    # Transform variables from radians to degrees
    for sym in [:pa, :ι, :χ, :η, :θo, :θs]
        if sym == :pa
            temp[:, indexin([sym,], prior_keys)[1]] .= (360 .+ temp[:, indexin([sym,], prior_keys)[1]] .* 180 / π  ) .% 360
        else
            temp[:, indexin([sym,], prior_keys)[1]] .= temp[:, indexin([sym,], prior_keys)[1]] .* 180 / π  
        end

    end
    Chains(temp, prior_keys)[[:m_d, :spin, :θo, :θs, :pa, :rpeak]]
end

include("pairplots_extensions.jl")
markersize = 2.5
theme_curr = Theme(
    Axis=(
        xticklabelsize=30,
        xlabelsize=50,
        yticklabelsize=30,
        ylabelsize=50,
        xticklabelrotation=pi / 4,
        xminorticks=IntervalsBetween(4),
        xlabelpadding=20,
        xticklabelfont="Computer Modern Roman",
        ylabelpadding=20,
        yticklabelfont="Computer Modern Roman",
    ),
    LineElement=(
        linewidth=10,
    ),
)

set_theme!(merge(theme_curr, theme_latexfonts()))
plt = Figure(resolution=(1150, 1150));
gs = GridLayout(plt[1:6, 1:6])
pairplot(
    gs,
    samples_to_plot => (
        PairPlots.Scatter(color=blue_cb, rasterize = true),#strokecolor = :black),
        MarginMakieHist(; bins=10, color=blue_cb, strokecolor=:black, rasterize = true),
        PairPlots.Contour(color=:black, bandwidth=2.0,rasterize = true),
    ),
    PairPlots.Truth(
        true_vals,
        color=:black,
        linewidth=8.0,
    ),
    PairPlots.Truth(
        nxcorr_vals[[:m_d, :spin, :θo, :θs, :rpeak]],
        color=red_cb,
        linewidth=3.0,
    ),
    axis=(;
        m_d=(;
            ticks=([2.0, 4.0, 6.0]),
            lims=(; low=prior.m_d.a, high=prior.m_d.b),
        ),
        spin=(;
            ticks=([-0.85, -0.45, -0.05]),
            lims=(; low=prior.spin.a, high=prior.spin.b),
        ),
        θo=(;
            lims=(; low=prior.θo.a*180/π, high=prior.θo.b*180/π),
            ticks=([10, 20, 30]),
        ),
        θs=(;
            lims=(; low=prior.θs.a*180/π, high=prior.θs.b*180/π),
            ticks=([50, 60, 70,80]),
        ),
        pa=(;
            lims=(; low=180, high=360),
            ticks=([200, 250, 300, 350]),
        ),
        rpeak=(;
            lims=(; low=prior.rpeak.a, high=prior.rpeak.b),
            ticks=([5, 10, 15]),
        ),),
    labels=Dict(:m_d => L"\theta_g", :spin => L"a", :θo => L"\theta_o", :pa => L"p.a.", :rpeak => L"R", :θs => L"\theta_s"),
);
label_models = [L"\textit{m}\text{F-ring}", L"\text{xs-ringauss}", L"\text{Hybrid Themage}"]
markers = [
    LineElement(color=:black, linestyle=:solid, linewidth=5),
    LineElement(color=red_cb, linestyle=:solid, linewidth=5),
]
labels = [L"\text{Truth}", L"\text{NxCORR Best Fit}"]

Legend(gs[4, 6],
    labelsize=33,
    [markers,],
    [labels,],
    [nothing,],
    tellheight=false,
    tellwidth=false,
    margin=(10, 10, 10, 10),
    halign=:right, valign=:bottom, orientation=:vertical,
)
colgap!(gs, 0)
rowgap!(gs, 0)
display(plt)

using Clustering
clusters = Clustering.kmeans(samples_to_plot.value[:,1:2,1]', 2)
cluster1 = samples_to_plot[clusters.assignments .== 1, :, :]
cluster2 = samples_to_plot[clusters.assignments .== 2, :, :]
fig = Figure();
ax = Axis(fig[1,1])
CairoMakie.scatter!(ax, cluster1.value[:,1:2,1])
CairoMakie.scatter!(ax, cluster2.value[:,1:2,1])
display(fig)
hpd(cluster1)
hpd(cluster2)

plt = Figure(resolution=(1150, 1150));
gs = GridLayout(plt[1:6, 1:6])
pairplot(
    gs,
    cluster1 => (
        #PairPlots.Scatter(color=blue_cb, rasterize = true),
        #PairPlots.HexBin(color=blue_cb, rasterize = true),
        #PairPlots.HexBin(colormap=CairoMakie.Makie.cgrad([white_cb_t, blue_cb]), rasterize=true),
        #PairPlots.Contour(color=:black, bandwidth=2.0, rasterize = true),
        PairPlots.Contourf(color=blue_cb_t1, bandwidth=3.0, rasterize = true),
        MarginMakieHist(; bins=10, color=blue_cb_t, strokecolor=:black, rasterize = true),
    ),
    cluster2 => (
        #PairPlots.Scatter(color=orange_cb, rasterize = true),
        #PairPlots.HexBin(colormap=CairoMakie.Makie.cgrad([white_cb_t, orange_cb]), rasterize=true),
        #PairPlots.Contour(color=:black, bandwidth=2.0, rasterize = true),
        PairPlots.Contourf(color=orange_cb_t1, bandwidth=3.0, rasterize = true),
        MarginMakieHist(; bins=10, color=orange_cb_t, strokecolor=:black, rasterize = true),
    ),
    PairPlots.Truth(
        true_vals,
        color=:black,
        linewidth=8.0,
    ),
    PairPlots.Truth(
        nxcorr_vals[[:m_d, :spin, :θo, :θs, :rpeak]],
        color=red_cb,
        linewidth=3.0,
    ),
    axis=(;
        m_d=(;
            ticks=([2.0, 4.0, 6.0]),
            lims=(; low=prior.m_d.a, high=prior.m_d.b),
        ),
        spin=(;
            ticks=([-0.85, -0.45, -0.05]),
            lims=(; low=prior.spin.a, high=prior.spin.b),
        ),
        θo=(;
            lims=(; low=prior.θo.a*180/π, high=prior.θo.b*180/π),
            ticks=([10, 20, 30]),
        ),
        θs=(;
            lims=(; low=prior.θs.a*180/π, high=prior.θs.b*180/π),
            ticks=([50, 60, 70,80]),
        ),
        pa=(;
            lims=(; low=180, high=360),
            ticks=([200, 250, 300, 350]),
        ),
        rpeak=(;
            lims=(; low=prior.rpeak.a, high=prior.rpeak.b),
            ticks=([5, 10, 15]),
        ),
        p1=(;
            lims=(; low=prior.p1.a, high=prior.p1.b),
            ticks=([2, 4, 6, 8]),
        ),
        p2=(;
            lims=(; low=prior.p2.a, high=prior.p2.b),
            ticks=([2, 4, 6, 8]),
        ),
        χ=(;
            lims=(; low=prior.χ.a*180/π, high=prior.χ.b*180/π),
        ),
        ι=(;
            lims=(; low=0*180/π, high=prior.ι.b*180/π),
            ticks=([0, 30, 60]),

        ),
        βv=(;
            lims=(; low=prior.βv.a, high=prior.βv.b),
            ticks=([0.1, 0.40, 0.8]),
        ),
        spec=(;
            lims=(; low=prior.spec.a, high=prior.spec.b),
            ticks=([1, 2, 3]),
        ),
        η=(;
            lims=(; low=prior.η.a*180/π, high=prior.η.b*180/π),
            ticks=([-150, -75, 0, 75, 150]),
        ),

    ),
    labels=Dict(:m_d => L"\theta_g", :spin => L"a", :θo => L"\theta_o", :pa => L"p.a.", :rpeak => L"R", :θs => L"\theta_s", :p1 => L"p_1", :p2 => L"p_2", :χ => L"χ", :ι => L"ι", :βv => L"β_v", :spec => L"\sigma", :η => L"η"), 
);
label_models = [L"\textit{m}\text{F-ring}", L"\text{xs-ringauss}", L"\text{Hybrid Themage}"]
markers = [
    LineElement(color=:black, linestyle=:solid, linewidth=5),
    LineElement(color=red_cb, linestyle=:solid, linewidth=5),
]
labels = [L"\text{Truth}", L"\text{NxCORR Best Fit}"]

Legend(gs[4, 6],
    labelsize=33,
    [markers,],
    [labels,],
    [nothing,],
    tellheight=false,
    tellwidth=false,
    margin=(10, 10, 10, 10),
    halign=:right, valign=:bottom, orientation=:vertical,
)
colgap!(gs, 0)
rowgap!(gs, 0)
display(plt)


save(joinpath((@__DIR__), "$(model_name)_reduced_posterior.pdf"), plt)
