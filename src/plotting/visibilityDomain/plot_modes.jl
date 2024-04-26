import CairoMakie as CMk
using Pyehtim
using Printf
using Comrade
using LaTeXStrings
using ImageFiltering

img1 = ehtim.image.load_fits("/n/home06/dochang/bamextension/visibilityDomain/plotting/mode1.fits").regrid_image(μas2rad(120), 180)
img1vals = reshape(pyconvert(Vector{Float64}, img1.imvec),(180,180))[180:-1:1,180:-1:1]

img2 = ehtim.image.load_fits("/n/home06/dochang/bamextension/visibilityDomain/plotting/mode2.fits").regrid_image(μas2rad(120), 180)
img2vals = reshape(pyconvert(Vector{Float64}, img2.imvec),(180,180))[180:-1:1,180:-1:1]

imgtrue = ehtim.image.load_fits("/n/home06/dochang/bamextension/visibilityDomain/plotting/truth.fits").regrid_image(μas2rad(120), 180)
imgtruevals = reshape(pyconvert(Vector{Float64}, imgtrue.imvec),(180,180))[180:-1:1,180:-1:1]

curr_theme = CMk.Theme(
    Axis = (
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
        xreversed = true,
        aspect = 18/14,
    ),
    Text = (
        fontsize = 25,
        font = "Computer Modern Roman"
    ),
    Heatmap = (
        colormap = :afmhot,
        colorbar = false,
        rasterize = true
    )
)
CMk.set_theme!(curr_theme)
fig = CMk.Figure(;figure_padding=5, resolution=((18/14*600)+10,400+10));
ax1 = CMk.Axis(fig[1,1]);
CMk.heatmap!(ax1, imgtruevals, colorrange=(0,max(imgtruevals...)*1.2))
CMk.text!(ax1, 170, 150, text=LaTeXString("Truth"), align=(:left, :top), color=:white)
CMk.text!(ax1, 170, 130, text = LaTeXString("\$θ_g=$(@sprintf("%.2f",3.83))\\;\\mu as\$"), align = (:left, :top), color = :white)
ax11 = CMk.Axis(fig[2,1]);
smimgtruevals = imfilter(imgtruevals, Kernel.gaussian((20*180/160)/(2*sqrt(2*log(2)))))
CMk.heatmap!(ax11, smimgtruevals, colorrange=(0,max(smimgtruevals...)*1.2))

ax2 = CMk.Axis(fig[1,2]);
CMk.heatmap!(ax2, img1vals, colorrange=(0,max(img1vals...)*1.2))
CMk.text!(ax2, 170, 150, text=LaTeXString("Cluster 1"), align=(:left, :top), color=:white)
CMk.text!(ax2, 170, 130, text = LaTeXString("\$θ_g=$(@sprintf("%.2f",3.629))\\;\\mu as\$"), align = (:left, :top), color = :white)
CMk.text!(ax2, 170, 45, text = LaTeXString("\$r_{peak}=$(@sprintf("%.f",3.75))\\;GM/c^2\$"), align = (:left, :top), color = :white)
ax12 = CMk.Axis(fig[2,2]);
smimg1vals = imfilter(img1vals, Kernel.gaussian((20*180/160)/(2*sqrt(2*log(2)))))
CMk.heatmap!(ax12, smimg1vals, colorrange=(0,max(smimg1vals...)*1.2))

ax3 = CMk.Axis(fig[1,3]);
CMk.heatmap!(ax3, img2vals, colorrange=(0,max(img2vals...)*1.2))
CMk.text!(ax3, 170, 150, text= LaTeXString("Cluster 2"), align=(:left, :top), color=:white)
CMk.text!(ax3, 170, 130, text = LaTeXString("\$θ_g=$(@sprintf("%.2f",1.6147))\\;\\mu as\$"), align = (:left, :top), color = :white)
CMk.text!(ax3, 170, 45, text = LaTeXString("\$r_{peak}=$(@sprintf("%.f",13.924))\\;GM/c^2\$"), align = (:left, :top), color = :white)

ax13 = CMk.Axis(fig[2,3]);
smimg2vals = imfilter(img2vals, Kernel.gaussian((20*180/160)/(2*sqrt(2*log(2)))))
CMk.heatmap!(ax13, smimg2vals, colorrange=(0,max(smimg2vals...)*1.2))

CMk.colgap!(fig.layout, 1, CMk.Fixed(0))
CMk.colgap!(fig.layout, 2, CMk.Fixed(0))
CMk.rowgap!(fig.layout, 1, CMk.Fixed(0))

CMk.ylims!(ax1, (20,160))
CMk.ylims!(ax2, (20,160))
CMk.ylims!(ax3, (20,160))
CMk.ylims!(ax11, (20,160))
CMk.ylims!(ax12, (20,160))
CMk.ylims!(ax13, (20,160))

display(fig)
CMk.save((@__DIR__)*"/modes.pdf", fig)

