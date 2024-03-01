using Revise
using Comrade
using Pyehtim
using Plots
using Pigeons
using LinearAlgebra;
using CairoMakie
using PairPlots
using MCMCChains
using ColorSchemes
using FileIO
using Krang
using Serialization
include(joinpath(dirname(@__DIR__) , "..", "models", "JuKeBOX.jl"))

red_cb = colorant"rgba(84%, 11%, 38%, 1.0)";
blue_cb = colorant"rgba(12%, 53%, 89%, 0.5)";
orange_cb = colorant"rgba(100%, 75%, 3%, 1.0)";
green_cb = colorant"rgba(0%, 30%, 25%, 1.0)"

model_name = "Data_2017"
in_base = abspath(dirname(@__DIR__), "..", "visibilityDomain","results","all","2024-01-31-23-56-33-MmVVnfWg")
include(joinpath(in_base, "model_params.jl"))
outpath = abspath(dirname(@__DIR__), ".." , "..","..","..","runs","visibilityDomain","data","$(model)")
#pt = PT(in_base)
#s_array = sample_array(pt)
rnd = 16
reduced_recorders = Serialization.deserialize(joinpath(in_base, "round=$rnd", "checkpoint", "reduced_recorders.jls"))
s_array = hcat([reduced_recorders.traces[n_tempering_levels=>i] for i in 1:2^rnd]...)'[:,1:end]

Plots.scatter(s_array[:, 1, :])
best_fit = begin
    param_file = open(abspath(dirname(@__DIR__), "..", "..","runs","image_domain","sa+0.94_r160_nall_tavg.fits", "JBOX", "best_nxcorr.txt"))
    [readline(param_file) for _ in 1:3]
    eval(Meta.parse(readline(param_file)[12:end]))
end

nxcorr_vals = (best_fit..., pa=360-72)
true_vals = (m_d=3.83, spin=-0.94, θo=17 / 180 * π, pa=360 - 72)

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

hpdi = hpd(samples_to_plot)

# the region to highlight will be a 50% quantile interval
intervals = collect(zip(hpdi.nt.lower, hpdi.nt.upper))

include(abspath((@__DIR__)*"/pairplots_extensions.jl"))
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
        linewidth=7,
    ),
)

set_theme!(merge(theme_curr, theme_latexfonts()))

modl = skymodel(post, transform(cpost, s_array[rand(1:size(s_array)[1]),:,1]))
modl |> Plots.plot
smoothed(modl, μas2rad(10)/(2√(2*log(2)))) |> Plots.plot
Plots.plot(dcphase)

modlphase = extract_table(obs, ClosurePhases())
modlphase.data.measurement .= closure_phases(modl, dcphase.config)  
modlamp = extract_table(obs, LogClosureAmplitudes())
modlamp.data.measurement .= logclosure_amplitudes(modl, dlcamp.config) 
chi2(modl, dcphase)/(length(dcphase))
chi2(modl, dlcamp)/(length(dlcamp))
Plots.plot(dcphase)
Plots.plot!(modlphase)
Plots.plot(dlcamp)
Plots.plot!(modlamp)

plt = Figure(resolution=(1150, 1150));
gs = GridLayout(plt[1:6, 1:6])
pairplot(gs,
    samples_to_plot => (
        PairPlots.Scatter(color=blue_cb, rasterize = true),
        PairPlots.Contour(color=:black, rasterize = true),
        MarginMakieHist(; bins=10, color=blue_cb, rasterize = true)
    ),
    PairPlots.Truth(
        (;
            m_d=2.05
        ),
        color=blue_cb,
        linewidth=4.5,
    ),
    PairPlots.Truth(
        (;
            m_d=3.62
        ),
        color=red_cb,
        linewidth=5.5,
    ),
    PairPlots.Truth(
        (;
            m_d=5.72
        ),
        color=green_cb,
        linewidth=5.5,
    ),
    PairPlots.Truth(
        (;
            m_d=3.83,
            pa = 360-72,
            θo=17 ,
        ),
        color=:black,
        linewidth=4.5,
        linestyle=:dash,
    ),
    PairPlots.Truth(
        (;
            m_d=3.16
        ),
        color=orange_cb,
        linewidth=5.5,
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
            ticks=([45,60,75]),
        ),
        pa=(;
            lims=(; low=180, high=360),
            ticks=([200, 250, 300, 350]),
        ),
        rpeak=(;
            lims=(; low=prior.rpeak.a, high=prior.rpeak.b),
            ticks=([5, 10, 15]),
        ),),
    labels=Dict(:m_d => L"\theta_g", :spin => L"a", :θo => L"\theta_o", :pa => L"p.a.", :rpeak => L"R", :θs => L"\theta_s", :p1 => L"p_1", :p2 => L"p_2", :χ => L"χ", :ι => L"ι", :βv => L"β_v", :spec_index => L"\sigma", :cross_spec_index => L"\sigma_{\zeta}", :η => L"η"),
)
colgap!(gs, 0)
rowgap!(gs, 0)
display(plt)
markers = [
    LineElement(color=blue_cb, linestyle=:solid, linewidth=5),
    LineElement(color=green_cb, linestyle=:solid, linewidth=5),
    LineElement(color=red_cb, linestyle=:solid, linewidth=5),
    LineElement(color=:black, linestyle=:dot, linewidth=5),
    LineElement(color=orange_cb, linestyle=:solid, linewidth=5),
]
labels = [L"\text{Walsh et al.}", L"\text{Simon et al.}", L"\text{Gebhardt et al.}", L"\text{EHT Collab.}", L"\text{Liepold et al.}"]

Legend(gs[4, 6],
    labelsize=40,
    [markers,],
    [labels,],
    [nothing,],
    tellheight=false,
    tellwidth=false,
    margin=(10, 10, 10, 10),
    halign=:right, valign=:bottom, orientation=:vertical,
)

display(plt)

posterior_sample = begin
    temp = tsamples.value[rand(1:size(tsamples)[1]), :, 1]
    NamedTuple{Tuple(temp.axes[1].val)}(temp.data)
end

posterior_sample = s_array[rand(1:size(s_array)[1]),:,1]
logdensityof(post, transform(cpost, posterior_sample))
logdensityof(post, rand(post.prior))
ax = Axis(gs[1:2, 3:4], aspect=1, yreversed=true)
hidedecorations!(ax)
april6 = ehtim.image.load_fits(joinpath((@__DIR__),"3598_blur_avg.fits")).regrid_image(μas2rad(135), 400)
april6.display()
april6 = reshape(pyconvert(Vector{Float64}, april6.imvec), (400,400))
CairoMakie.heatmap!(
    ax, 
    april6, 
    colormap=:afmhot,
    colorrange=(0, 3.50e-5),
    rasterize=true
    )
text!(ax, 10, 45; text="April 6, 2017", fontsize=30, color=:white, depth_shift=-1.0f0)

smoothedmdlimage = smoothed(skymodel(post, transform(cpost, posterior_sample)), μas2rad(20 / (2√(2log(2)))))
smoothedintmap = intensitymap(smoothedmdlimage, μas2rad(165), μas2rad(165), 400, 400) |> Array
ax = Axis(gs[1:2, 5:6], aspect=1, xreversed=true)
hidedecorations!(ax)
CairoMakie.heatmap!(
    ax,
    smoothedintmap,
    colormap=:afmhot,
    colorrange=(0, 5.35e-4), 
    rasterize=true
)
CairoMakie.xlims!(ax, 350, 50)
CairoMakie.ylims!(ax, 28, 328)
text!(ax, 340, 290; text=L"20\;μ\text{as blurred sample}", fontsize=30, color=:white)
display(plt)

save((@__DIR__) * "/$(model_name)_reduced_posterior.pdf", plt)

