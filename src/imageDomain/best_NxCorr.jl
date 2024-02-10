using VIDA
using Plots
using StatsBase
using Comrade
using OptimizationBBO
using ImageFiltering

include(joinpath((@__DIR__),"..", "models", "JuKeBOX.jl"))
include(joinpath((@__DIR__),"..","models","GpuizedModel.jl"))
include(joinpath((@__DIR__),"..","models","GpuizedBlurredModel.jl"))
include(joinpath((@__DIR__),"..","models","defaults.jl"))
include(abspath(joinpath((@__DIR__),"..", "models","vidawrappers.jl")))
inbase = abspath((@__DIR__),"..","..", "data", "GRMHD")
readdir(inbase)

maxevals=100_000
blur  = 10f0
dual_cone_blur_temp = create_dual_cone_blur_model(blur)
files = String[] 
for data_file in readdir(inbase)
    if occursin("M_a", data_file)
        file_name = joinpath(inbase, data_file)
        append!(files, [file_name])
    end
end
model_info_list = ModelInfo.(Ref(dual_cone_blur_temp), Ref("JBOXBlur"), Ref(lower_dc), Ref(upper_dc), files, blur)
fit.(model_info_list, maxevals)