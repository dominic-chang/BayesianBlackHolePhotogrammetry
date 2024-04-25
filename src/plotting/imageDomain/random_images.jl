using Enzyme
using CairoMakie
using Comrade
using Plots
using Krang
using Zygote
Enzyme.Compiler.RunAttributor[] = false
Enzyme.Compiler.CheckNan[] = false
include(joinpath(dirname(@__DIR__) , "..", "models", "JuKeBOX.jl"))

red_cb = colorant"rgba(84%, 11%, 38%, 1.0)";
blue_cb = colorant"rgba(12%, 53%, 89%, 0.5)";
orange_cb = colorant"rgba(100%, 75%, 3%, 1.0)";
green_cb = colorant"rgba(0%, 30%, 25%, 1.0)"


lines = readlines(abspath(joinpath(dirname(@__DIR__), ".." , "..", "runs", "image_domain","M_a-0.94_Rh160_i30.fits","JBOX","best_nxcorr.txt")))
eval(Meta.parse(lines[end]))
function custom_dualcone(θ, metadata)
    #(;nmax, cache) = metadata
    (;nmax, ) = metadata
    params = (
        spin = θ.spin,
        θo   = θ.f*θ.θs,
        θs   = θ.θs,
        rpeak= θ.rpeak,
        p1   = θ.p1,
        p2   = θ.p1,
        χ    = θ.χ,
        ι    = θ.ι,
        βv   = θ.βv,
        spec    = θ.spec,
        η    = θ.η
    )
    model = JuKeBOX(nmax, params)
    return stretched(model, μas2rad(4.0), μas2rad(4.0))
    return Comrade.modelimage(rotated(stretched(model, μas2rad(θ.m_d), μas2rad(θ.m_d)), θ.pa), cache, true)
end
function dualcone(θ, metadata)
    #(;nmax, cache) = metadata
    (;nmax, ) = metadata
    model = JuKeBOX(nmax, θ)
    return stretched(model, μas2rad(4.0), μas2rad(4.0))
    return Comrade.modelimage(rotated(stretched(model, μas2rad(θ.m_d), μas2rad(θ.m_d)), θ.pa), cache, true)
end
metadat = (nmax = 2, cache=Comrade.create_cache(FFTAlg(), IntensityMap(zeros(100,100), μas2rad(120), μas2rad(120))))

using Distributions, VLBIImagePriors, VLBISkyModels
prior = (;
    #m_d = Uniform(4.0, 4.01),
    spin = Uniform(-1,-0.9999),
    θo =Uniform(30/180*π, 31/180*π),
    #f = Uniform(0.01, 1.0),
    θs =Uniform(1/180*π, 90/180*π),
    pa = Uniform(-π, π),
    rpeak= Uniform(1., 10.),
    p1 = Uniform(0.1, 10),
    p2 = Uniform(2, 10),
    χ = Uniform(-π, π),
    ι = Uniform(0.0, π/2),
    βv = Uniform(0.0 ,0.99),
    spec = Uniform(-1.,3.),
    η = Uniform(0,π),
)
ndims = length(prior)
pprior = VLBIImagePriors.NamedDist(prior)
cprior = ascube(pprior)
grid = imagepixels(μas2rad(120), μas2rad(120), 150, 150)


using CairoMakie
curr_theme = Theme(
    Axis=(
        aspect=1,
        xticksvisible=false, 
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        spinewidth=4
        ),
    Heatmap = (
        colormap=:afmhot,
        rasterize=false
    )
)
CairoMakie.set_theme!(merge(curr_theme, theme_latexfonts()))
fig = Figure(size= (800,420));
gs = CairoMakie.GridLayout(fig[1:4, 1:8])
axs = [Axis(gs[i,j]) for i in 1:4 for j in 1:8]
for ax in axs
    params = transform(cprior, rand(ndims))
    if params.θo > params.θs
        ax.topspinecolor= red_cb
        ax.bottomspinecolor= red_cb
        ax.leftspinecolor= red_cb
        ax.rightspinecolor= red_cb
    else
        ax.topspinecolor= blue_cb
        ax.bottomspinecolor= blue_cb
        ax.leftspinecolor= blue_cb
        ax.rightspinecolor= blue_cb
    end
    m = smoothed(dualcone(params, metadat), μas2rad(2))

    collect(intensitymap(m, grid))
    CairoMakie.heatmap!(ax, collect(intensitymap(m, grid)))
end
CairoMakie.colgap!(gs, 4)
CairoMakie.rowgap!(gs, 0)
display(fig)
save((@__DIR__) * "/random.pdf", fig)
