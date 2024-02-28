preamble_path = joinpath((@__DIR__),"dual_cone_data_2017_preamble.jl")
include(preamble_path)

using Pigeons 

settings = Pigeons.MPISettings(;
submission_system=:slurm, 
add_to_submission = [
    "#SBATCH -p wholenode",
    ], 
    environment_modules=["intel/19.0.5.281","impi/2019.5.281"]
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
    	walltime="01:00:00",
        n_threads = 12,
        dependencies = [
            Pigeons, # <- Pigeons itself can be skipped, added automatically
	    preamble_path
        ]
    ),
    n_rounds=14
)
