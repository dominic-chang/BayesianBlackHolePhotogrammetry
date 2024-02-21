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
blur = 0f0#μas2rad(10.0f0)
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

model_info_list = ModelInfo.(Ref(eq_blur_temp), Ref("EqBlur"), Ref(lower_eq), Ref(upper_eq), files, blur)
(; model, name, lower, upper, file, blur) = model_info_list[1]
inimg = (VIDA.load_image(file))
dx, dy = Float32.(([pixelsizes(inimg)...]))
g = axisdims(inimg)
vals = zeros(Float32, size(inimg))

if iszero(blur)
    vals = Float32.(parent(inimg))
else
    σ_px = blur / (2 * sqrt(2 * log(2f0))) / abs(dx)
    σ_py = blur / (2 * sqrt(2 * log(2f0))) / abs(dy)

    # Now I need to pick my kernel size. I am going out to 5σ for the
    # gaussian kernel. I have to add one for the convolution to play nice
    nkern = Int(floor(σ_px) * 10 + 1)
    vals = Float32.(imfilter(parent(inimg),
        Kernel.gaussian((σ_py, σ_px), (nkern, nkern)),
        Fill(0.0, inimg),
        Algorithm.FFT()
    ))
end
g32 = RectiGrid((X=map(Float32, (g.X)), Y=map(Float32, (g.Y))))
img = IntensityMap(vals, g32) #|>MtlArray
Plots.plot(img)

bh = VIDA.NxCorr(img)

best_fit = (m_d = (5.017628061882903), spin = -0.9346459134142151, θo = 0.4785851172566213, θs = 1.263293226080634, rpeak = 4.9218205590076565, p1 = 0.10032446760976753, p2 = 3.7916779586163902, χ = -3.1415788543946332, ι = 1.5707637112185893, βv = 0.17565582217360834, σ = 0.0003075599473727954, η = -0.05925432839814393, x0 = (0.4432871400486391), y0 = (0.013455392422866197))
mdlimg = intensitymap(dual_cone_blur_temp(best_fit), imagepixels(μas2rad(200f0), μas2rad(200f0), 400, 400))
mdlimg  |> Plots.plot
nxcorr(img, mdlimg)
divergence(bh, dual_cone_blur_temp(best_fit))
best_fit = (m_d = 5.340221491299042, spin = -0.9997162334137694, θo = 0.3104509369404855, θs = 0.992095997114991, rpeak = 6.278945282186813, p1 = 0.10001432508035828, p2 = 4.040149068160272, χ = 0.10742861074813215, ι = 1.2182960980284654, βv = 0.4236918288930219, σ = 1.0263133775893358, η = 0.13422924352941967, x0 = 1.9989785678084822, y0 = 0.3319156921225219)
mdlimg = intensitymap(dual_cone_blur_temp(best_fit), imagepixels(μas2rad(200f0), μas2rad(200f0), 400, 400))
mdlimg  |> Plots.plot
nxcorr(img, mdlimg)
using Cthulhu
divergence(bh, dual_cone_blur_temp(best_fit))
bh.mimg |> Plots.plot
divergence(bh, dual_cone_temp(best_fit))
#intensitymap!(bh.mimg, dual_cone_blur_temp(best_fit)) |> Plots.plot
#intensitymap!(bh.mimg, dual_cone_temp(best_fit)) |> Plots.plot
bh.mimg |> Plots.plot
divergence(bh, dual_cone_temp(best_fit))
best_fit = (m_d = 5.226105463961742, spin = -0.9906787902909712, θo = 0.3137702455121262, θs = 0.9791688826230474, rpeak = 6.776765321260125, p1 = 0.13456024590004964, p2 = 4.193689888445819, χ = 0.041275693956445814, ι = 1.1805560517044458, βv = 0.4066724700370679, σ = 0.5932285608045533, η = 0.15109669244258894, x0 = 1.9561785751642322, y0 = 0.45214669041690847)
mdlimg  |> Plots.plot
nxcorr(img, mdlimg)
divergence(bh, dual_cone_blur_temp(best_fit))
model_info_list = ModelInfo.(Ref(eq_temp), Ref("Eq"), Ref(lower_eq), Ref(upper_eq), files, 0.0f0)
imgfit.(model_info_list, maxevals)

#model_info_list = ModelInfo.(Ref(eq_blur_temp), Ref("EqBlur"), Ref(lower_eq), Ref(upper_eq), files, blur)
#imgfit.(model_info_list, maxevals)

model_info_list = ModelInfo.(Ref(dual_cone_temp), Ref("JBOX"), Ref(lower_dc), Ref(upper_dc), files, 0.0f0)
imgfit.(model_info_list, maxevals)

#model_info_list = ModelInfo.(Ref(dual_cone_blur_temp), Ref("JBOXBlur"), Ref(lower_dc), Ref(upper_dc), files, blur)
#imgfit.(model_info_list, maxevals)

model_info_list = ModelInfo.(Ref(eq_dual_cone_temp), Ref("EqDualCone"), Ref(lower_eq_dc), Ref(upper_eq_dc), files, 0.0f0)
imgfit.(model_info_list, maxevals)

#model_info_list = ModelInfo.(Ref(eq_dual_cone_blur_temp), Ref("EqDualConeBlur"), Ref(lower_eq_dc), Ref(upper_eq_dc), files, blur)
#imgfit.(model_info_list, maxevals)


