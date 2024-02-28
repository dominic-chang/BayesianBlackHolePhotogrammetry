using Comrade
using Pyehtim
using Plots
using Pigeons
using LinearAlgebra; LinearAlgebra.BLAS.set_num_threads(1)
using FFTW; FFTW.set_num_threads(1)
using Krang
include((@__DIR__)*"/../models/JuKeBOX.jl")

sze = 180
n_tempering_levels = 20
modelfov = 120
seed = 4
fractional_noise = 0.01
phasecal = true
ampcal = true
add_th_noise = true
scan_avg = true
nmax = 1

function dualcone(θ, metadata)
    (;nmax, cache) = metadata
    model = JuKeBOX(nmax, θ)
    return Comrade.modelimage(rotated(stretched(model, μas2rad(θ.m_d), μas2rad(θ.m_d)), θ.pa), cache, true)
end

obs = ehtim.obsdata.load_uvfits(joinpath((@__DIR__),"..","..","data","2017","SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
obs = scan_average(obs.flag_uvdist(uv_min=0.1e9))
obs = obs.add_fractional_noise(0.01)
dlcamp = extract_table(obs, LogClosureAmplitudes())
dcphase = extract_table(obs, ClosurePhases())

using Distributions
prior = (;
    m_d = Uniform(1.5, 8.0),
    spin = Uniform(-1,-0.01),
    θo =Uniform(1/180*π, 40/180*π),
    θs =Uniform(40/180*π, 90/180*π),
    pa = Uniform(-π, π),
    rpeak= Uniform(1., 10.),
    p1 = Uniform(0.1, 10),
    p2 = Uniform(1, 10),
    χ = Uniform(-π, π),
    ι = Uniform(0.0, π/2),
    βv = Uniform(0.0 ,0.99),
    spec = Uniform(-1.,3.),
    η = Uniform(0,π),
)
cache = create_cache(NFFTAlg(dlcamp), IntensityMap(zeros(sze, sze), (μas2rad(modelfov)), (μas2rad(modelfov))))
metadat = (nmax=nmax, cache=cache)
lklhd = RadioLikelihood(dualcone, dlcamp, dcphase ;skymeta=metadat)

using VLBIImagePriors
ndim = length(prior) # gets diminsions of parameter space
pprior = VLBIImagePriors.NamedDist(prior)
post = Posterior(lklhd, pprior)
cpost = ascube(post)
log_posterior = Comrade.logdensityof(cpost)
cprior= Distributions.product_distribution(fill(Uniform(), ndim))
log_prior = Pigeons.DistributionLogPotential(cprior)

cprsample = rand(cprior)
prsample = transform(cpost, cprsample)
model = dualcone(prsample, metadat)
model |> Plots.plot