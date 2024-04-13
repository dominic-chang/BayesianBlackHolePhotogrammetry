preamble_path = joinpath((@__DIR__),"dual_cone_data_2017_preamble.jl")
include(preamble_path)

using Pigeons 

append!(Pigeons._rosetta.custom,["srun",`sbatch`,`scancel`,"#SBATCH", "--job-name=","-o ", "-e ", "\$SLURM_SUBMIT_DIR", `squeue --job`, `squeue -u`, `sinfo`])
function Pigeons.resource_string(m::MPIProcesses, ::Val{:custom})
    return """
    #SBATCH -t $(m.walltime)
    #SBATCH --ntasks=$(m.n_mpi_processes)
    #SBATCH --cpus-per-task=$(m.n_threads)
    #SBATCH --mem-per-cpu=$(m.memory)
    """
end

settings = Pigeons.MPISettings(;
submission_system=:slurm, 
add_to_submission = [
    "#SBATCH -p blackhole",
    ], 
    environment_modules=["intel","intelmpi"]
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
        ],
        mpiexec_args=`--mpi=pmi2`
    ),
    n_rounds=14
)
