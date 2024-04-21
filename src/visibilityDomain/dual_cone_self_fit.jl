preamble_path = joinpath((@__DIR__),"dual_cone_self_fit_preamble.jl")

include(preamble_path)

using Pigeons 

params= (
    exec = "srun",
    submit = `sbatch`,
    del = `scancel`,
    directive = "#SBATCH",
    job_name = "--job-name=",
    output_file = "-o ",
    error_file = "-e ",
    submit_dir = "\$SLURM_SUBMIT_DIR",
    job_status = `squeue --job`,
    job_status_all = `squeue -u`,
    ncpu_info = `sinfo`
)

Pigeons.add_custom_submission_system(params)

function Pigeons.resource_string(m::Pigeons.MPIProcesses, ::Val{:custom})
    return """
    #SBATCH -t $(m.walltime)
    #SBATCH --ntasks=$(m.n_mpi_processes)
    #SBATCH --cpus-per-task=$(m.n_threads)
    #SBATCH --mem=32G
    """
end

settings = Pigeons.MPISettings(;
submission_system=:custom, 
add_to_submission = [
    "#SBATCH -p blackhole",
    ], 
    environment_modules=["intel","intelmpi"]
)
Pigeons.setup_mpi(settings)

#pt = PT(abspath(dirname(@__DIR__), ".." ,"results", "latest"))
#pt = increment_n_rounds!(pt, 4)
#
#result = pigeons(pt.exec_folder, 
#    on = Pigeons.MPIProcesses(
#        n_mpi_processes = n_tempering_levels,
#    	walltime="10-00:00:00",
#        n_threads = 48,
#        dependencies = [
#            Pigeons, # <- Pigeons itself can be skipped, added automatically
#	    preamble_path
#        ],
#        mpiexec_args=`--mpi=pmi2`
#    )
#)


pt = Pigeons.pigeons(
    target=log_posterior, 
    reference=log_prior; 
    record = [traces, round_trip, Pigeons.timing_extrema], 
    checkpoint=true, 
    n_chains=n_tempering_levels, 
    on = Pigeons.MPIProcesses(
        n_mpi_processes = n_tempering_levels,
    	walltime="14-00:00:00",
        n_threads = 24,
        dependencies = [
            Pigeons, # <- Pigeons itself can be skipped, added automatically
	    preamble_path
        ],
        mpiexec_args=`--mpi=pmi2`
    ),
    n_rounds=20
)
