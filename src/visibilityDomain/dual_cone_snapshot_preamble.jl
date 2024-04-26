using Comrade
using Pyehtim
using Plots
using Pigeons
using LinearAlgebra; LinearAlgebra.BLAS.set_num_threads(1)
using FFTW; FFTW.set_num_threads(1)
using Krang
include(joinpath(dirname(@__DIR__), "models", "JuKeBOX.jl"))

sze = 50
n_tempering_levels = 22
modelfov = 100
seed = 4
fractional_noise = 0.01
rotate_fits = 108/180*π
phasecal = true
ampcal = true
add_th_noise = true
scan_avg = true
nmax = 1

function dualcone(θ, metadata)
    (;nmax, cache) = metadata
    model = JuKeBOX(nmax, θ)
    mdl = Comrade.modelimage(rotated(stretched(model, μas2rad(θ.m_d), μas2rad(θ.m_d)), θ.pa), cache, true)
    return mdl /flux(mdl)
end
inbase = abspath(dirname(@__DIR__), "..", "data")
obsin = ehtim.obsdata.load_uvfits(joinpath(inbase, "2017", "SR1_M87_2017_096_lo_hops_netcal_StokesI.uvfits"))
inimg = ehtim.image.load_image(joinpath(inbase, "GRMHD", "snapshots", "image_ma+0.5_1000_163_20_nall.h5"))
inimg.rf = obsin.rf
inimg.ra = obsin.ra
inimg.dec = obsin.dec
inimg = inimg.rotate(rotate_fits)
obs = inimg.observe_same(obsin, ampcal=ampcal, phasecal=phasecal, add_th_noise=add_th_noise, seed=seed, ttype="fast")
obs = scan_average(obs.flag_uvdist(uv_min=0.1e9))
obs = obs.add_fractional_noise(fractional_noise)

dlcamp = extract_table(obs, LogClosureAmplitudes())
dcphase = extract_table(obs, ClosurePhases())

using Distributions
prior = (;
    m_d = Uniform(1.5, 8.0),
    spin = Uniform(-1,-0.01),
    θo =Uniform(1/180*π, 40/180*π),
    θs =Uniform(40/180*π, 90/180*π),
    pa = Uniform(-π, 0),
    rpeak= Uniform(1., 18.),
    p1 = Uniform(0.1, 10),
    p2 = Uniform(1, 10),
    χ = Uniform(-π, π),
    ι = Uniform(0.0, π/2),
    βv = Uniform(0.0 ,0.99),
    spec = Uniform(-1.,3.),
    η = Uniform(-π,π),
)
cache = create_cache(NFFTAlg(dlcamp), IntensityMap(zeros(sze, sze), (μas2rad(modelfov)), (μas2rad(modelfov))))
metadat = (nmax=nmax, cache=cache)
lklhd = RadioLikelihood(dualcone, dlcamp, dcphase ;skymeta=metadat)

using VLBIImagePriors
ndim = length(prior) # gets diminsions of parameter space
pprior = VLBIImagePriors.NamedDist(prior)
post = Posterior(lklhd, pprior)
cpost = ascube(post)
log_posterior = cpost
cprior= Distributions.product_distribution(fill(Uniform(), ndim))
log_prior = Pigeons.DistributionLogPotential(cprior)

cprsample = rand(cprior)
prsample = transform(cpost, cprsample)
model = dualcone(prsample, metadat)
#model |> Plots.plot
