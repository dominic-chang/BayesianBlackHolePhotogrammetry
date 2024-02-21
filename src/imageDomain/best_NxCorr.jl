using Revise 
using VIDA
using Plots
using StatsBase
using Comrade
using OptimizationBBO
using ImageFiltering
using Krang
using VLBIImagePriors
using ComradeBase

include(joinpath((@__DIR__), "..", "models", "Equatorial.jl"))
include(joinpath((@__DIR__), "..", "models", "JuKeBOX.jl"))
include(joinpath((@__DIR__), "..", "models", "EqDualCone.jl"))
include(joinpath((@__DIR__), "..", "models", "GpuizedModel.jl"))
include(joinpath((@__DIR__), "..", "models", "GpuizedBlurredModel.jl"))
include(joinpath((@__DIR__), "..", "models", "defaults.jl"))
include(abspath(joinpath((@__DIR__), "..", "models", "vidawrappers.jl")))
inbase = abspath((@__DIR__), "..", "..", "data", "GRMHD")
readdir(inbase)

maxevals = 100_000
blur = μas2rad(10.0f0)
eq_blur_temp = create_eq_blur_model(blur)
dual_cone_blur_temp = create_dual_cone_blur_model(blur)
eq_dual_cone_blur_temp = create_eq_dual_cone_blur_model(blur)

files = String[]
for data_file in readdir(inbase)
    if occursin("M_a", data_file)
        file_name = joinpath(inbase, data_file)
        append!(files, [file_name])
    end
end

model_info_list = ModelInfo.(Ref(eq_temp), Ref("Eq"), Ref(lower_eq), Ref(upper_eq), files, 0.0f0)
imgfit.(model_info_list, maxevals)

model_info_list = ModelInfo.(Ref(eq_blur_temp), Ref("EqBlur"), Ref(lower_eq), Ref(upper_eq), files, blur)
imgfit.(model_info_list, maxevals)

model_info_list = ModelInfo.(Ref(dual_cone_temp), Ref("JBOX"), Ref(lower_dc), Ref(upper_dc), files, 0.0f0)
imgfit.(model_info_list, maxevals)

model_info_list = ModelInfo.(Ref(dual_cone_blur_temp), Ref("JBOXBlur"), Ref(lower_dc), Ref(upper_dc), files, blur)
imgfit.(model_info_list, maxevals)

model_info_list = ModelInfo.(Ref(eq_dual_cone_temp), Ref("EqDualCone"), Ref(lower_eq_dc), Ref(upper_eq_dc), files, 0.0f0)
imgfit.(model_info_list, maxevals)

model_info_list = ModelInfo.(Ref(eq_dual_cone_blur_temp), Ref("EqDualConeBlur"), Ref(lower_eq_dc), Ref(upper_eq_dc), files, blur)
imgfit.(model_info_list, maxevals)


