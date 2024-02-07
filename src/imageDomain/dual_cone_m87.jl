include((@__DIR__)*"/../models/JuKeBOX.jl")
include((@__DIR__)*"/../models/GpuizedModel.jl")
using VIDA

function JuKeBOX(m_d::T,spin::T,θo::T,θs::T,rpeak::T,p1::T,p2::T,χ::T,ι::T,βv::T,spec::T,η::T, x0::T,y0::T) where T
    return modify(JuKeBOX(spin,θo,θs,rpeak,p1,p2,χ,ι,βv,spec, η, 2),
    Stretch(m_d, m_d), Shift(x0, y0))
end

using VLBISkyModels
function dualCone_temp(θ)
    return GPUModel(
        JuKeBOX(
        θ.m_d,
        θ.spin,
        θ.θo,
        θ.θs,
        θ.rpeak,
        θ.p1,
        θ.p2,
        θ.χ,
        θ.ι,
        θ.βv,
        θ.σ,
        θ.η,
        θ.x0,
        θ.y0
    )
    )
end

using Plots
img = VIDA.load_image((@__DIR__)*"/../../data/GRMHD/M_a-0.94_Rh160_i30.fits")
Plots.plot(img)
struct NxCORR{T} <: VIDA.AbstractDivergence
    img::T
    mimg::T
end

bh =  VIDA.NxCorr(img)
using StatsBase

start = (
    m_d =μas2rad(4.920864394023698e0),
    spin = -0.9342827323237523e0,
    θo = 0.4711974848553482e0,
    θs = 1.3518693322434683e0,
    rpeak = 5.59261856012756e0,
    p1 = 0.158968560099307e0,
    p2 = 4.192568509151882e0,
    χ = 2.7242920822576653e0,
    ι = 1.0113763707982746e0,
    βv = 0.3096525555378556e0,
    σ = -0.34430247937749187e0,
    η = 0.01817958773540323e0,
    x0 = μas2rad(0.7511569255298416e0),
    y0 = μas2rad(0.18446651454923932e0)
 )

dualCone_temp(start)|> x->intensitymap(x, μas2rad(120e0), μas2rad(120e0), 200, 200) |> Plots.plot
img = dualCone_temp(start)|> x->intensitymap(x, μas2rad(120e0), μas2rad(120e0), 200, 200) 
img = img ./ sum(img)
sum(img)/length(img)
divergence(bh, dualCone_temp(start))
lower = (
    m_d = μas2rad(2.0e0),
    spin = -1.0e0,
    θo = 0e0/180*π,
    θs = π/4e0,
    rpeak = 1e0,
    p1 = 0.1e0,
    p2 = 1.0e0,
    χ = -1.0e0π,
    ι = 0.0e0,
    βv = 0.01e0,
    σ= 0.0e0,
    η = -1.0e0π,
    x0 = -μas2rad(10e0),
    y0 = -μas2rad(10e0)
)
dualCone_temp(lower)|> x->intensitymap(x, μas2rad(120e0), μas2rad(120e0), 400, 400) |> Plots.plot
divergence(bh, dualCone_temp(lower))
upper = (
    m_d = μas2rad(8.0),
    spin = -0.01,
    θo = 89.9/180*π,
    θs = π/2,
    rpeak = 18.0,
    p1 = 10.0,
    p2 = 10.0,
    χ =1.0π,
    ι = π/2,
    βv = 0.999,
    σ= 10.0,
    η = 1.0π,
    x0 = μas2rad(10.0),
    y0 = μas2rad(10.0)
)
dualCone_temp(upper) |> x->intensitymap(x, μas2rad(120), μas2rad(120), 40, 40) |> Plots.plot
divergence(bh, dualCone_temp(upper))
prob = VIDAProblem(bh, dualCone_temp, lower, upper);

using OptimizationBBO

keys(upper)
count = 0::Int
f, t, (lb, ub) = VIDA.build_opt(prob, true)
callback = function (state, loss; doplot = false)
    global count
    display(count)
    count += 1
    return false
end
@time xopt, optfilt, divmin = vida(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=50_000, callback=callback)
triptic(img, optfilt)


#using Enzyme
#template = dualCone_temp(start)
#test_func(x,y) = let model = template
#    return VIDA.intensity_point(model, (X=x, Y=y))
#end
#
#test_func(μas2rad(15.0), μas2rad(5.0))
#
#devs = []
#αvals = -100:1.1:100
#βvals = -100:1.1:100
#
#α =rand(αvals)
#β = rand(βvals)
#Enzyme.autodiff(Enzyme.Reverse, test_func, Active, Active(μas2rad(α)), Active(μas2rad(β)))