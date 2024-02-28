include((@__DIR__)*"/dual_cone_data_2017_preamble.jl")

using Pigeons 

settings = Pigeons.MPISettings(;
submission_system=:slurm, 
add_to_submission = [
    "#SBATCH -p blackhole",
    ], 
    environment_modules=["intel/23.2.0-fasrc01","intelmpi/2021.10.0-fasrc01"]
)
Pigeons.setup_mpi(settings)

pt = Pigeons.pigeons(
    target=log_posterior, 
    reference=log_prior; 
    record = [traces, round_trip, Pigeons.timing_extrema], 
    checkpoint=true, 
    n_chains=n_tempering_levels, 
    on = Pigeons.MPIProcesses(
        n_mpi_processes = n_tempering_levels,
        n_threads = 48,
        dependencies = [
            Pigeons, # <- Pigeons itself can be skipped, added automatically
            (@__DIR__)*"/DualCone_Data_2017_preamble.jl"
        ]
    ),
    n_rounds=18
)
