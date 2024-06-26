import CairoMakie as CMk
using Comrade
using PythonCall
using Meshes
import Krang
using Colors
using ColorSchemes
using Pyehtim
using Printf
using LaTeXStrings
include(joinpath((@__DIR__), "..", "..", "models", "customModels.jl"))

red_cb = CMk.colorant"rgba(84%, 11%, 38%, 1.0)";
blue_cb = CMk.colorant"rgba(12%, 53%, 89%, 1.0)";
orange_cb = CMk.colorant"rgba(100%, 75%, 3%, 1.0)";

curr_theme = CMk.Theme(
    Axis=(
        xtickcolor=:white,
        xtickalign=1,
        xticksmirrored=true,
        xticklabelsvisible=false,
        xgridvisible=false,
        ytickcolor=:white,
        ytickalign=1,
        yticksmirrored=true,
        yticklabelsvisible=false,
        ygridvisible=false,
        leftspinecolor=:white,
        rightspinecolor=:white,
        topspinecolor=:white,
        bottomspinecolor=:white,
    ),
    Lines=(
        linewidth=5,
    ), 
    Text=(
        fontsize=20,
    )
)
CMk.set_theme!(merge(curr_theme, CMk.theme_latexfonts()))

c1 = colorant"rgba(0%, 0%, 0%, 1.0)";
c2 = orange_cb;
c3 = blue_cb;
c4 = red_cb;
oranges = range(c1, stop=c2, length=150);
blues = range(c1, stop=c3, length=150);
reds = range(c1, stop=c4, length=150);

models = ["ma+0.5_r20_nall_tavg","sa+0.5_r160_nall_tavg", "ma+0.94_r10_nall_tavg","sa+0.94_r160_nall_tavg"]

model = models[1]
inimg = ehtim.image.load_fits(joinpath((@__DIR__) , "..","..","..","data","GRMHD", "$model.fits"))
pyconvert(Vector{Float64}, inimg.imvec)
np = pyimport("numpy")

function doublePower(rs, p1, p2)
    return r -> (r/rs)^(p1) / (1 + (r/rs)^(p1+p2))
end

fig = CMk.Figure(resolution=(995, 660), figure_padding=1);
for (i, model) in enumerate(models)
    np = pyimport("numpy")

    inimg = ehtim.image.load_fits(joinpath((@__DIR__) , "..","..","..","data","GRMHD", "$model.fits"))
    inimg = inimg.rotate((180 + 108) / 180 * π)
    inimg = inimg.regrid_image(μas2rad(120), 500)
    ax1 = CMk.Axis(fig[1, i], aspect=1.5, xreversed=true)
    CMk.ylims!(ax1, 83, 416)
    CMk.heatmap!(ax1, reshape(pyconvert(Vector{Float64}, inimg.imvec), (500, 500)), colormap=:afmhot, rasterize=true)
    tm = CMk.text!(ax1, 485, 400, text=CMk.L"\text{GRMHD (truth)}", align=(:left, :top), color=:white)
    if i == 1
        CMk.linesegments!(ax1, [20, 187], [150, 150], color=:white)
        CMk.text!(ax1, 485, 360, text=CMk.L"\text{MAD}", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 320, text=CMk.L"R_{\text{high}}=20", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 270, text=CMk.L"\theta_o=17^\circ", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 220, text=CMk.L"a=0.5", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 170, text=CMk.L"\theta_g=3.83\;\mu\text{as}", align=(:left, :top), color=:white)
        CMk.text!(ax1, 135, 140, text=CMk.L"40\,\mu as", align=(:left, :top), color=:white)
    elseif i == 2
        CMk.text!(ax1, 485, 360, text=CMk.L"\text{SANE}", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 320, text=CMk.L"R_{\text{high}}=160", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 270, text=CMk.L"\theta_o=17^\circ", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 220, text=CMk.L"a=0.5", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 170, text=CMk.L"\theta_g=3.83 \mu\text{as}", align=(:left, :top), color=:white)
    elseif i == 3
        CMk.linesegments!(ax1, [20, 187], [150, 150], color=:white)
        CMk.text!(ax1, 485, 360, text=CMk.L"\text{MAD}", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 320, text=CMk.L"R_{\text{high}}=10", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 270, text=CMk.L"\theta_o=17^\circ", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 220, text=CMk.L"a=0.94", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 170, text=CMk.L"\theta_g=3.83\;\mu\text{as}", align=(:left, :top), color=:white)
        CMk.text!(ax1, 135, 140, text=CMk.L"40\,\mu as", align=(:left, :top), color=:white)
    elseif i == 4
        CMk.text!(ax1, 485, 360, text=CMk.L"\text{SANE}", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 320, text=CMk.L"R_{\text{high}}=160", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 270, text=CMk.L"\theta_o=17^\circ", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 220, text=CMk.L"a=0.94", align=(:left, :top), color=:white)
        CMk.text!(ax1, 485, 170, text=CMk.L"\theta_g=3.83 \mu\text{as}", align=(:left, :top), color=:white)

    end
    

    # Plot best fit curve EqModel
    lines = readlines(joinpath(dirname(@__DIR__), ".." , "..", "runs", "image_domain","$model.fits","Eq","best_nxcorr.txt"))
    eval(Meta.parse(lines[end]))
    inimg = ehtim.image.load_fits(joinpath(dirname(@__DIR__), ".." , "..","runs","image_domain","$model.fits","Eq","best.fits"))
    inimg = inimg.rotate((180 + 108) / 180 * π)
    inimg = inimg.regrid_image(μas2rad(120), 500)
    ax2 = CMk.Axis(fig[2, i], aspect=1.5, xreversed=true, spinewidth=2)
    CMk.ylims!(ax2, 83, 416)
    CMk.heatmap!(ax2, reshape(pyconvert(Vector{Float64}, inimg.imvec), (500, 500)), colormap=:afmhot, rasterize=true)
    CMk.text!(ax2, 485, 380, text=CMk.L"\text{Equatorial}", align=(:left, :top), color=:white)
    CMk.text!(ax2, 485, 330, text = LaTeXString("\$θ_s=90^\\circ\$"), align = (:left, :top), color = :white)
    CMk.text!(ax2, 485, 280, text = LaTeXString("\$θ_o=$(@sprintf("%.f",best_fit.θo*180/π))^\\circ\$"), align = (:left, :top), color = :white)
    CMk.text!(ax2, 485, 230, text = LaTeXString("\$a=$(@sprintf("%.2f",-best_fit.spin))"), align = (:left, :top), color = :white)
    CMk.text!(ax2, 485, 180, text = LaTeXString("\$θ_g=$(@sprintf("%.2f",best_fit.m_d)) \\;\\mu\$as"), align = (:left, :top), color = :white)
    θs = π / 2
    rs = best_fit.rpeak
    p1 = best_fit.p1
    p2 = best_fit.p2
    prof = doublePower(rs, p1, p2)

    xvals = sin(θs) .* (1.25:30)
    yvals = cos(θs) .* (1.25:30)
    colorvals = prof.(sqrt.(xvals .^ 2 .+ yvals .^ 2))

    # Plot best fit curve DualCone
    lines = readlines(abspath(joinpath(dirname(@__DIR__), ".." , "..", "runs", "image_domain","$model.fits","JBOX","best_nxcorr.txt")))
    eval(Meta.parse(lines[end]))
    inimg = ehtim.image.load_fits(joinpath(dirname(@__DIR__), ".." , "..","runs","image_domain","$model.fits","JBOX","best.fits"))
    inimg = inimg.rotate((180 + 108) / 180 * π)
    inimg = inimg.regrid_image(μas2rad(120), 500)
    ax3 = CMk.Axis(fig[3, i], aspect=1.5, xreversed=true, spinewidth=2)
    CMk.ylims!(ax3, 83, 416)
    CMk.heatmap!(ax3, reshape(pyconvert(Vector{Float64}, inimg.imvec), (500, 500)), colormap=:afmhot, rasterize=true)
    CMk.text!(ax3, 485, 380, text=CMk.L"\text{DualCone}", align=(:left, :top), color=:white)
    CMk.text!(ax3, 485, 330, text = LaTeXString("\$θ_s=$(@sprintf("%.f",best_fit.θs*180/π))^\\circ\$"), align = (:left, :top), color = :white)
    CMk.text!(ax3, 485, 280, text = LaTeXString("\$θ_o=$(@sprintf("%.f",best_fit.θo*180/π))^\\circ\$"), align = (:left, :top), color = :white)
    CMk.text!(ax3, 485, 230, text = LaTeXString("\$a=$(@sprintf("%.2f",-best_fit.spin))"), align = (:left, :top), color = :white)
    CMk.text!(ax3, 485, 180, text = LaTeXString("\$θ_g=$(@sprintf("%.2f",best_fit.m_d)) \\;\\mu\$as"), align = (:left, :top), color = :white)
    θs = best_fit.θs
    rs = best_fit.rpeak
    p1 = best_fit.p1
    p2 = best_fit.p2
    prof = doublePower(rs, p1, p2)

    xvals = sin(θs) .* (1.25:30)
    yvals = cos(θs) .* (1.25:30)
    colorvals = prof.(sqrt.(xvals .^ 2 .+ yvals .^ 2))

    # Plot best fit curve EqDualCone
    lines = readlines(abspath(joinpath(dirname(@__DIR__), ".." , "..", "runs", "image_domain","$model.fits","EqDualCone","best_nxcorr.txt")))
    eval(Meta.parse(lines[end]))
    inimg = ehtim.image.load_fits(joinpath(dirname(@__DIR__), ".." , "..","runs","image_domain","$model.fits","EqDualCone","best.fits"))
    inimg = inimg.rotate((180 + 108) / 180 * π)
    inimg = inimg.regrid_image(μas2rad(120), 500)
    ax4 = CMk.Axis(fig[4, i], aspect=1.5, xreversed=true, spinewidth=2)
    CMk.ylims!(ax4, 83, 416)
    CMk.heatmap!(ax4, reshape(pyconvert(Vector{Float64}, inimg.imvec), (500, 500)), colormap=:afmhot, rasterize=true)
    CMk.text!(ax4, 485, 380, text=CMk.L"\text{Equatorial + Dual Cone}", align=(:left, :top), color=:white)
    CMk.text!(ax4, 485, 330, text = LaTeXString("\$θ_s=$(@sprintf("%.f",best_fit.θs*180/π))^\\circ\$"), align = (:left, :top), color = :white)
    CMk.text!(ax4, 485, 280, text = LaTeXString("\$θ_o=$(@sprintf("%.f",best_fit.θo*180/π))^\\circ\$"), align = (:left, :top), color = :white)
    CMk.text!(ax4, 485, 230, text = LaTeXString("\$a=$(@sprintf("%.2f",-best_fit.spin))"), align = (:left, :top), color = :white)
    CMk.text!(ax4, 485, 180, text = LaTeXString("\$θ_g=$(@sprintf("%.2f",best_fit.m_d)) \\;\\mu\$as"), align = (:left, :top), color = :white)
    θs = best_fit.θs
    rs = best_fit.r_cone
    p1 = best_fit.p1_cone
    p2 = best_fit.p2_cone
    prof = doublePower(rs, p1, p2)

    xvals = sin(θs) .* (1.25:20)
    yvals = cos(θs) .* (1.25:20)
    colorvals = prof.(sqrt.(xvals .^ 2 .+ yvals .^ 2))

end
CMk.colgap!(fig.layout, CMk.Fixed(0.1))
CMk.rowgap!(fig.layout, CMk.Fixed(0.1))
display(fig)

CMk.save((@__DIR__) * "/m87_nxcorr_summary.pdf", fig)
