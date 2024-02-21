using Revise 
using VIDA
using CairoMakie
using StatsBase
using Comrade
using OptimizationBBO
using ImageFiltering
using Printf
using Pyehtim
using Krang
using VLBIImagePriors
using ComradeBase
using LaTeXStrings
include(joinpath((@__DIR__), "..", "..","models", "JuKeBOX.jl"))
include(joinpath((@__DIR__), "..", "..","models", "Equatorial.jl"))
include(joinpath((@__DIR__), "..", "..","models", "EqDualCone.jl"))
include(joinpath((@__DIR__), "..", "..", "models", "GpuizedModel.jl"))
include(joinpath((@__DIR__), "..", "..", "models", "GpuizedBlurredModel.jl"))
include(joinpath((@__DIR__), "..", "..", "models", "defaults.jl"))
include(abspath(joinpath((@__DIR__), "..", "..", "models", "vidawrappers.jl")))
inbase = abspath((@__DIR__), "..", "..", "data", "GRMHD")
red_cb = colorant"rgba(84%, 11%, 38%, 1.0)";
blue_cb = colorant"rgba(12%, 53%, 89%, 1.0)";
orange_cb = colorant"rgba(100%, 75%, 3%, 1.0)";
green_cb = colorant"rgba(0%, 3%, 25%, 1.0)";

curr_theme = Theme(
	Axis = (
		#xgridvisible=false,
		#xticksvisible=false,
		xticklabelsvisible = false,
		xtickcolor = :white,
		xtickalign = 1,
		xticksmirrored = true,
		#ygridvisible=false,
		#yticksvisible=false,
        #xticklabelfont="New Computer Modern",
        #yticklabelfont="New Computer Modern",
        #ylabelfont="New Computer Modern",
		yticklabelsvisible = false,
		ytickcolor = :white,
		ytickalign = 1,
		yticksmirrored = true,
	),
)
sze = 400
set_theme!(merge(curr_theme, theme_latexfonts()))
inbase = abspath(joinpath((@__DIR__) , "..", "..", "..", "data", "GRMHD"))
dirlist = filter(x -> x[1] == 'M', readdir(abspath(joinpath((@__DIR__) , "..", "..", "..", "runs", "image_domain"))))

fig = CairoMakie.Figure(resolution = (1000*0.8, 650*0.8), figure_padding = 2);
for (count, currdirname) in enumerate(dirlist)

	ax1 = CairoMakie.Axis(fig[1, count], aspect = 2)
	ax2 = CairoMakie.Axis(fig[2, count], aspect = 2)
	ax3 = CairoMakie.Axis(fig[3, count], aspect = 2)
	ax4 = CairoMakie.Axis(fig[4, count], aspect = 2)

	ylims!(ax1, 100, 300); ylims!(ax2, 100, 300); ylims!(ax3, 100, 300);ylims!(ax4, 100, 300)


	currimg = Pyehtim.ehtim.image.load_fits(joinpath(inbase , currdirname ))

	inimage = VIDA.load_image(abspath(joinpath(inbase , currdirname )))
	img = reverse(reshape(pyconvert(Vector{Float64}, currimg.ivec), (400, 400)), dims = 2)
	hm = CairoMakie.heatmap!(ax1, img, colormap = :afmhot, rasterize=true)
	text!(ax1, 10, 290, text = LaTeXString("\$θ_o=$(currdirname[17:18])^\\circ\$"), align = (:left, :top), color = :white)
	text!(ax1, 10, 250, text = LaTeXString("\$θ_g=5.03"), align = (:left, :top), color = :white)
	text!(ax1, 10, 210, text = LaTeXString("\$a=0.94"), align = (:left, :top), color = :white)

	plot_best_fit(param_file, model, ax) = let
		curr_file = open(param_file, "r")
		[readline(curr_file) for _ in 1:3]
		ans = eval(Meta.parse(readline(curr_file)[12:end]))
		close(curr_file)
		best_fit = NamedTuple{(propertynames(ans))}([Float32(i) for i in ans])

		start = model(best_fit)
		fit_img = intensitymap(start, imagepixels(pyconvert.(Ref(Float32), (currimg.fovx(), currimg.fovy()))...,sze,sze))
		CairoMakie.heatmap!(ax, reverse(Array(fit_img), dims=1), colormap = :afmhot, rasterize=true)
		text!(ax, 10, 290, text = LaTeXString("\$θ_o=$(@sprintf("%.f",best_fit.θo*180/π))^\\circ\$"), align = (:left, :top), color = :white)
		text!(ax, 10, 250, text = LaTeXString("\$θ_g=$(@sprintf("%.2f",best_fit.m_d))"), align = (:left, :top), color = :white)
		text!(ax, 10, 210, text = LaTeXString("\$a=$(@sprintf("%.2f",-best_fit.spin))"), align = (:left, :top), color = :white)
	end


	param_file1 = abspath(joinpath((@__DIR__) , "..", "..", "..", "runs", "image_domain" , currdirname , "Eq", "best_nxcorr.txt"))
	best_fit = plot_best_fit(param_file1, eq_temp, ax2)
	param_file2 = abspath(joinpath((@__DIR__) , "..", "..", "..", "runs", "image_domain" , currdirname , "JBOX", "best_nxcorr.txt"))
	best_fit = plot_best_fit(param_file2, dual_cone_temp, ax3)
	param_file3 = abspath(joinpath((@__DIR__) , "..", "..", "..", "runs", "image_domain" , currdirname , "EqDualCone", "best_nxcorr.txt"))
	best_fit = plot_best_fit(param_file3, eq_dual_cone_temp, ax4)

	if count == 1

		ax1.ylabel = "GRMHD\n(truth)"
		ax2.ylabel = "Equatorial"
		ax3.ylabel = "Dual Cone"
		ax4.ylabel = "Equatorial\n+ Dual Cone"
		linesegments!(ax4, [30, 110], [145, 145], color = :white)
		tm = text!(ax4, 30, 140, text = L"40\,\mu as", align = (:left, :top), color = :white)
	end
end

indir = joinpath((@__DIR__) , "..", "..","..","runs", "image_domain") |> abspath
directory_names = filter(x -> occursin("M_a", x), readdir(indir))
inclinations = map(x -> x[end-1:end], directory_names)
models = ["JBOX", "Eq", "EqDualCone"]#, "EqDualConeBlur"]

nxcorr_vals = Dict([model => Vector{Float64}() for model in models])

for (inc, dir) in zip(inclinations, directory_names)
	for model in models
		f = open(joinpath(indir, dir, model, "best_nxcorr.txt"))

		append!(nxcorr_vals[model], [parse(Float64, readline(f)[begin+9:end])])
		close(f)
	end
end

ax = Axis(
	fig[5, 1:4],
	ylabel = LaTeXString("NxCORR"),
	xlabel = LaTeXString("Observer Inclination (\${}^\\circ\$)"),
	xticklabelsvisible = true,
	xticks = ([1, 2, 3, 4], ["30", "50", "70", "90"]),
	yticks = [0.85, 0.9, 0.95, 1.0],
	yticklabelsvisible = true,
	xlabelsize=20
);
xlims!(ax, (0.5, 4.5))
ylims!(ax, (0.84, 1.005))

colors = [red_cb, blue_cb, orange_cb]#, :blue]
elements = [:circle, :rect, :utriangle]
model_names = reverse(collect(keys(nxcorr_vals)))
for (color, model, element) in zip(colors, model_names, elements)
	scatter!(ax, nxcorr_vals[model], color = color, marker = element)
end
model_markers = [MarkerElement(marker = element, color = color,
	strokecolor = :transparent,
	markersize = 20) for (color, element) in zip(colors, elements)]
model_labels = [L"\text{Equatorial}",L"\text{Equatorial + Dual Cone}", L"\text{Dual Cone}", ]
Legend(fig[5, 1:4], model_markers, model_labels, valign = :bottom, halign = :left)

colgap!(fig.layout, 3)
rowgap!(fig.layout, 1)
display(fig)

save(joinpath((@__DIR__), "nxcorr_summary.pdf"), fig)
