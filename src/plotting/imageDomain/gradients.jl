using Enzyme
using CairoMakie
using Comrade
using Plots
using Krang
using Zygote
Enzyme.Compiler.RunAttributor[] = false
Enzyme.Compiler.CheckNan[] = false
include(joinpath(dirname(@__DIR__) , "..", "models", "JuKeBOX.jl"))

lines = readlines(abspath(joinpath(dirname(@__DIR__), ".." , "..", "runs", "image_domain","M_a-0.94_Rh160_i30.fits","JBOX","best_nxcorr.txt")))
eval(Meta.parse(lines[end]))
function dualcone(θ, metadata)
    #(;nmax, cache) = metadata
    (;nmax, ) = metadata
    model = JuKeBOX(nmax, θ)
    return stretched(model, μas2rad(θ.m_d), μas2rad(θ.m_d))
    return Comrade.modelimage(rotated(stretched(model, μas2rad(θ.m_d), μas2rad(θ.m_d)), 0.0), cache, true)
end
metadat = (nmax = 2, cache=Comrade.create_cache(FFTAlg(), IntensityMap(zeros(100,100), μas2rad(120), μas2rad(120))))
best_fit = NamedTuple{map(x-> x == :σ ? :spec : x, keys(best_fit))}(values(best_fit))
m = dualcone(best_fit, metadat)

grid = imagepixels(μas2rad(120), μas2rad(120), 100, 100)
Plots.plot(intensitymap(m, grid))
img = intensitymap(m, grid) 
img |> Plots.plot
rimg = zero(img)
dimg = zero(img)
function point(m_d, spin, θo, θs, rpeak, p1, p2, χ, ι, βv, spec, η, X)
    best_fit = (m_d=m_d, spin=spin, θo=θo, θs=θs, rpeak=rpeak, p1=p1, p2=p2, χ=χ, ι=ι, βv=βv, spec=spec, η=η)
    model = dualcone(best_fit, (nmax=1,))
    return ComradeBase.intensity_point(model, (X=X[1], Y=X[2]))
end

npix = 400
fovxy = 120

vcat([
    collect(abs.(Enzyme.autodiff(Enzyme.Reverse, point, Active, Active.(values(best_fit)[begin:end-2])..., Const(μas2rad.([x,y])))[1][begin:end-1]))
    for x in range(-fovxy,fovxy, length=npix)
    for y in range(-fovxy/2,fovxy/2, length=npix)
]...)
grads = reshape(vcat([
    collect(abs.(Enzyme.autodiff(Enzyme.Reverse, point, Active, Active.(values(best_fit)[begin:end-2])..., Const(μas2rad.([x,y])))[1][begin:end-1]))
    for x in range(-fovxy,fovxy, length=npix)
    for y in range(-fovxy/2, fovxy/2, length=npix)
]...), (12, npix, npix))

curr_theme = CairoMakie.Theme(
    Axis=(
        xticksvisible=false,
        yticksvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        titlesize=25,
        xreversed=true,
        aspect=2,
        topspinevisible=false,
        bottomspinevisible=false,
        leftspinevisible=false,
        rightspinevisible=false,

    ),
    Text=(
        fontsize=30,
        textcolor=:white,
        color=:white
    ),
    Colorbar=(
        ticklabelsize=20,
        ticksvisible=false,
        flipaxis=true,
        vertical = false,
        width=400
        #tickwidth=2,
        #tickcolor=:white,
        #ticklabelcolor=:white,
        #titlecolor=:white,
        #titlefontsize=40,
        #titlecolor=:white
    ),
    Heatmap=(
        colormap=:afmhot,
        rasterize=true
    )
)
CairoMakie.set_theme!(merge(curr_theme, CairoMakie.theme_latexfonts()))
begin
fig = Figure(size=(1300,730))
gs = GridLayout(fig[1:3,1:5])
labels = reshape(collect(string.(keys(best_fit)[begin:end-2])), (4,3))
axs = reshape([Axis(gs[i, 1+j]) for i in 1:3 for j in 1:4], (3,4))

ax = Axis(fig[1:3, 1], aspect=2/3)
grid = imagepixels(μas2rad(fovxy), μas2rad(3fovxy/2), 400, 400)
img = intensitymap(m, grid) 
CairoMakie.heatmap!(ax, range(-fovxy,fovxy, length=npix), range(-fovxy,fovxy, length=npix), collect(img), xlims=(-fovxy,fovxy))
CairoMakie.text!(ax, -100, -110, text="Dual Cone", align=(:right,:bottom), justification=:right)
CairoMakie.arrows!(ax, [0,], [0,], [0,], [70,], color=:white, linewidth=3, arrowsize=20)


label_names = (m_d=L"m_d", spin=L"a", θo=L"θ_o", θs=L"θ_s", rpeak=L"R", p1=L"p_1", p2=L"p_2", χ=L"\chi", ι=L"\iota", βv=L"β_v", spec=L"\sigma", η=L"\eta")
for I in LinearIndices(axs)
    CairoMakie.heatmap!(axs[I], range(-fovxy,fovxy, length=npix), range(-fovxy/2,fovxy/2, length=npix), reverse(grads[I,:,:]',dims=1), xlims=(-fovxy,fovxy), colormap=:viridis)
    lstring = "\$\\left|\\frac{\\partial P_{ij}}{\\partial "*label_names[keys(best_fit)[I]].s[2:end-1]*"}\\right|\$"
    CairoMakie.text!(axs[I], -120, -50, text=CairoMakie.Makie.LaTeXString("$lstring"), align=(:right,:bottom), justification=:right, fontsize=22)
end
Colorbar(fig[4, 2:5], limits = (0, 10), colormap = :viridis, ticks=(0:5:10, ["min","", "max"]))

axlab = Axis(fig[5,1], topspinevisible=false, bottomspinevisible=false, leftspinevisible=false, rightspinevisible=false)
hidedecorations!(axlab)
CairoMakie.xlims!(axlab, 0, 1)
CairoMakie.ylims!(axlab, 0, 1)
CairoMakie.text!(axlab, 0.5, 0.5, text="Image", color=:black, justificaton=:center, align= (:center, :center))

axlab = Axis(fig[5,2:5], aspect=10)
hidedecorations!(axlab)
CairoMakie.xlims!(axlab, 0, 1)
CairoMakie.ylims!(axlab, 0, 1)
CairoMakie.text!(axlab, 0.5, 0.5, text=L"\text{Magnitude of Image Derivatives }\left(|\partial P_{ij}/\partial\theta|\right)",color=:black, justificaton=:center, align= (:center, :center))

lax = Axis(fig[0:5, 1:2], aspect=0.1, xticksvisible=false, yticksvisible=false, xticklabelsvisible=false, yticklabelsvisible=false)
hidedecorations!(lax)
CairoMakie.xlims!(lax, 0, 1)
CairoMakie.ylims!(lax, 0, 1)

CairoMakie.linesegments!(lax, [0.485, 0.485], [0, 1], color=:black, linewidth=10)  

rowgap!(gs, 2)
colgap!(gs, 2)
display(fig)
save((@__DIR__) * "/gradients.pdf", fig)


end
