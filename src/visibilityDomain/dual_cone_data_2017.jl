include((@__DIR__)*"/dual_cone_data_2017_preamble.jl")

using Pigeons 

pt = Pigeons.pigeons(
    target=log_posterior, 
    reference=log_prior; 
    record = [traces, round_trip, Pigeons.timing_extrema], 
    checkpoint=true, 
    n_chains=n_tempering_levels, 
    n_rounds=16
)
using MCMCChains
using CairoMakie
using PairPlots

samples = Chains(sample_array(pt), variable_names(pt))
transform(cpost, vec(samples.value[1,:,:]))
outchains = map(x->collect(transform(cpost,x)), [vec(samples.value[i,:,:]) for i in 1:1024])
tchains = reshape(transpose(hcat(outchains...)), (4096,24,1))
samples.value
#pairplot(keys(prior), samples.value)
pairplot(Chains(tchains, [i for i in keys(prior)]))
#sample(post,)
out_chains = Chains(tchains, [i for i in keys(prior)])
# Save a chain.
outpath = abspath((@__DIR__)*"/../runs/visibilityDomain/data6/m872017")
mkpath(outpath)
using Serialization
serialize(outpath*"/chain-file.jls", out_chains)

key_names = Tuple([i for i in keys(prior)])
peak_posterior = NamedTuple{key_names}(out_chains.value[end,:,:])
model = eqdualcone(peak_posterior, metadata)
Plots.plot(model)#comrade.intensitymap(zeros(sze, sze), g)

size(out_chains.value)
rand((1:size(out_chains.value)[1]))
rand_posterior = NamedTuple{key_names}(out_chains.value[rand((1:size(out_chains.value)[1])),:,:])
model = eqdualcone(rand_posterior, metadata)
Plots.plot(smoothed(model, μas2rad(10)))#comrade.intensitymap(zeros(sze, sze), g)


Plots.plot(model)#comrade.intensitymap(zeros(sze, sze), g)
Comrade.save(model,outpath*"/temp.fits")
