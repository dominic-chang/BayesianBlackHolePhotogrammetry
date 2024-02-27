import VIDA
using Plots
import CairoMakie as CMk
using StatsBase
using Comrade
using OptimizationBBO
using ImageFiltering
using Krang
using VLBIImagePriors
using ComradeBase
using Colors

blur = μas2rad(10f0)

include(joinpath((@__DIR__), "..", "..", "models", "customModels.jl"))
red_cb = colorant"rgba(84%, 11%, 38%, 1.0)";
blue_cb = colorant"rgba(12%, 53%, 89%, 1.0)";
orange_cb = colorant"rgba(100%, 75%, 3%, 1.0)";
green_cb = colorant"rgba(0%, 3%, 25%, 1.0)";


data_path = joinpath((@__DIR__), "..", "..", "..", "data", "GRMHD")
grmhd_model = "M_a-0.94_Rh160_i30"
inimg_path = joinpath(data_path, grmhd_model*".fits")
inimg = (VIDA.load_image(inimg_path))
dx, dy = Float32.(([pixelsizes(inimg)...]))
g = axisdims(inimg)
vals = zeros(Float32, size(inimg))

σ_px = blur / (2 * sqrt(2 * log(2f0))) / abs(dx)
σ_py = blur / (2 * sqrt(2 * log(2f0))) / abs(dy)

nkern = Int(floor(σ_px) * 10 + 1)
vals = Float32.(imfilter(parent(inimg),
    Kernel.gaussian((σ_py, σ_px), (nkern, nkern)),
    Fill(0.0, inimg),
    Algorithm.FFT()
))
g32 = RectiGrid((X=map(Float32, (g.X)), Y=map(Float32, (g.Y))))
img = IntensityMap(Float32.(parent(inimg)), g32)
img |> Plots.plot
img_blur = IntensityMap(vals, g32)
img_blur |> Plots.plot

results_path = joinpath(replace(replace(inimg_path, "data"=> "runs"), "GRMHD" => "image_domain"), "JBOX", "best_nxcorr.txt")
best_fit = begin
    infile = open(results_path, "r")
    for _ in 1:3
        readline(infile)
    end
    eval(Meta.parse(readline(infile)))
end

results_blur_path = replace(results_path, "JBOX"=> "JBOXBlur")
best_blur_fit = begin
    infile = open(results_blur_path, "r")
    for _ in 1:3
        readline(infile)
    end
    result = eval(Meta.parse(readline(infile)))
    NamedTuple{keys(result)}(Float32.(values(result)))
end

dual_cone_blur_temp = create_dual_cone_blur_model(blur)

model = dual_cone_temp_threaded(best_fit)
model_blur = dual_cone_blur_temp(best_blur_fit)
intmap_model = intensitymap(model, g)
intmap_model_blur = intensitymap(model_blur, g32)
resid =  100 .* ((img./flux(img)) .- (intmap_model./flux(intmap_model))) ./ max((intmap_model./flux(intmap_model))...)
resid_blur = 100 .* ((img_blur./flux(img_blur)) .- (intmap_model_blur./flux(intmap_model_blur))) ./ max((intmap_model_blur./flux(intmap_model_blur))...)

curr_theme = CMk.Theme(
    Axis=(
        aspect=1,
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        xgridvisible=false,
        ygridvisible=false,
        xticklabelsize=17,
        xlabelsize=20,
        ylabelsize=20
        ),
    Text=(
        fontsize=18,
        color=:white,
        ),
    Colorbar=(
        ticklabelsize=20,
        labelsize=20,
        ),
    Heatmap=(
        rasterize=true,
        ),
)
CMk.set_theme!(merge(curr_theme, CMk.theme_latexfonts()))
rbcustom = append!(
    collect(range(start=red_cb, stop=colorant"rgba(256,256,256,1.0)", length=19))[begin:end-1],
    range(start=colorant"rgba(256,256,256,1.0)", stop=blue_cb, length=19)
);

begin
    fig = CMk.Figure(size=(1200,670));
    ax1 = CMk.Axis(fig[1,1]);

    CMk.heatmap!(ax1, reverse(Array(img), dims=1), colormap=:afmhot, clims=(-50,50))
    dx, dy = pixelsizes(inimg)
    CMk.lines!(ax1, [10,10+μas2rad(40) / dx] , [50,50], color=:white, strokewidth=2)    
    CMk.text!(ax1, 10, 50, text="40 μas", align=(:left, :top), color=:white, fontsize=20)
    CMk.lines!(ax1, [0,400] , [200,200], color=blue_cb, linewidth=1.5)    
    CMk.lines!(ax1 , [200,200], [0,400], color=red_cb, linewidth=1.5)    
    CMk.text!(ax1, 390, 390, text="Truth", align=(:right, :top), color=:white)


    ax2 = CMk.Axis(fig[1,2]);
    CMk.heatmap!(ax2, reverse(Array(intmap_model), dims=1), colormap=:afmhot, clims=(-50,50))
    dx, dy = pixelsizes(intmap_model)
    CMk.lines!(ax2, [10,10+μas2rad(40) / dx] , [50,50], color=:white, strokewidth=2)    
    CMk.text!(ax2, 10, 50, text="40 μas", align=(:left, :top), color=:white, fontsize=20)
    CMk.lines!(ax2, [0,400] , [200,200], color=blue_cb, linewidth=3, linestyle=:dash)    
    CMk.lines!(ax2 , [200,200], [0,400], color=red_cb, linewidth=3, linestyle=:dash)    
    CMk.text!(ax2, 390, 390, text="Dual Cone", align=(:right, :top), color=:white)

    ax3 = CMk.Axis(fig[1,3], topspinevisible=false, rightspinevisible=false, xlabel="RA. DEC chords (μas)", ylabel="Intensity", xticksvisible=true, xticklabelsvisible=true);
    CMk.lines!(ax3, rad2μas.(collect(img.X)), reverse(img[:,200])./flux(img), color=blue_cb, linewidth=1.5)
    CMk.lines!(ax3, rad2μas.(collect(intmap_model.X)), reverse(intmap_model[:,200])./flux(intmap_model), color=blue_cb, linestyle=CMk.Linestyle([0.0, 2, 4, 5.0]), linewidth=3.0)
    CMk.lines!(ax3, rad2μas.(collect(img.Y)), reverse(img[200, :]./flux(img)), color=red_cb, linewidth=1.5)
    CMk.lines!(ax3, rad2μas.(collect(intmap_model.Y)), reverse(intmap_model[200, :]./flux(intmap_model)), color=red_cb, linestyle=CMk.Linestyle([0.0, 2, 4, 5.0]), linewidth=3.0)

    ax4 = CMk.Axis(fig[1,4]);
    hm = CMk.heatmap!(ax4, reverse(Array(resid), dims=1), colormap=reverse(rbcustom), clims=(-50,50), colorrange=(-30,30))
    CMk.Colorbar(fig[1,5], hm, label="Residual Difference (%)")    

    ax1 = CMk.Axis(fig[2,1]);

    CMk.heatmap!(ax1, reverse(Array(img_blur), dims=1), colormap=:afmhot, clims=(-50,50))
    dx, dy = pixelsizes(img_blur)
    CMk.lines!(ax1, [10,10+μas2rad(40) / dx] , [50,50], color=:white, strokewidth=2)    
    CMk.text!(ax1, 10, 50, text="40 μas", align=(:left, :top), color=:white, fontsize=20)
    lin_img_hor = CMk.lines!(ax1, [0,400] , [200,200], color=blue_cb, linewidth=1.5)    
    lin_img_ver = CMk.lines!(ax1 , [200,200], [0,400], color=red_cb, linewidth=1.5)    
    CMk.text!(ax1, 10, 390, text="10 μas Blur", align=(:left, :top), color=:white)
    CMk.text!(ax1, 390, 390, text="Truth", align=(:right, :top), color=:white)

    ax2 = CMk.Axis(fig[2,2]);
    CMk.heatmap!(ax2, reverse(Array(intmap_model_blur), dims=1), colormap=:afmhot, clims=(-50,50))
    dx, dy = pixelsizes(intmap_model_blur)
    CMk.lines!(ax2, [10,10+μas2rad(40) / dx] , [50,50], color=:white, strokewidth=2)    
    CMk.text!(ax2, 10, 50, text="40 μas", align=(:left, :top), color=:white, fontsize=20)
    lin_mdl_hor = CMk.lines!(ax2, [0,400] , [200,200], color=blue_cb, linewidth=3, linestyle=:dash)    
    lin_mdl_ver = CMk.lines!(ax2 , [200,200], [0,400], color=red_cb, linewidth=3, linestyle=:dash)    
    CMk.text!(ax2, 10, 390, text="10 μas Blur", align=(:left, :top), color=:white)
    CMk.text!(ax2, 390, 390, text="Dual Cone", align=(:right, :top), color=:white)

    ax3 = CMk.Axis(fig[2,3], topspinevisible=false, rightspinevisible=false, xlabel="RA. DEC chords (μas)", ylabel="Intensity", xticksvisible=true, xticklabelsvisible=true);
    CMk.lines!(ax3, rad2μas.(collect(img_blur.X)), reverse(img_blur[:,200])./flux(img_blur), color=blue_cb, linewidth=1.5)
    CMk.lines!(ax3, rad2μas.(collect(intmap_model_blur.X)), reverse(intmap_model_blur[:,200])./flux(intmap_model_blur), color=blue_cb, linestyle=CMk.Linestyle([0.0, 2, 4, 5.0]), linewidth=3.0)
    CMk.lines!(ax3, rad2μas.(collect(img_blur.Y)), reverse(img_blur[200, :]./flux(img_blur)), color=red_cb, linewidth=1.5)
    CMk.lines!(ax3, rad2μas.(collect(intmap_model_blur.Y)), reverse(intmap_model_blur[200, :]./flux(intmap_model_blur)), color=red_cb, linestyle=CMk.Linestyle([0.0, 2, 4, 5.0]), linewidth=3.0)

    ax4 = CMk.Axis(fig[2,4]);
    hm = CMk.heatmap!(ax4, reverse(Array(resid_blur), dims=1), colormap=reverse(rbcustom), clims=(-50,50), colorrange=(-6,6))
    CMk.Colorbar(fig[2,5], hm, label="Residual Difference (%)")    
    model_markers = [lin_img_hor, lin_img_ver, lin_mdl_hor, lin_mdl_ver]
    model_labels = ["Image Horizontal Section", "Image Vertical Section", "Model Horizontal Section", "Model Vertical Section"]
    CMk.Legend(fig[3, 1:5], model_markers, model_labels, valign = :top, halign = :center, labelsize=20, height=30, orientation=:horizontal)
    CMk.colgap!(fig.layout, 1)
    CMk.rowgap!(fig.layout, 1)


    display(fig)
    CMk.save(joinpath((@__DIR__), "triptic.pdf"),fig)

end