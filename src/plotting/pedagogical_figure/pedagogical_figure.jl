#using Pkg;
#Pkg.activate((@__DIR__));
import GLMakie as GLMk
using GeometryBasics, LinearAlgebra, Colors, ColorSchemes, Random, Meshes, Rotations, ImageCore, Krang, CoordinateTransformations
using FileIO
import Meshes
using Rotations
GLMk.Makie.inline!(true)

spin = 0.99

curr_theme = GLMk.Theme(
    Axis=(
        xgridvisible=false,
        ygridvisible=false,
        xspinesvisible=false,
        yspinesvisible=false,
        yticklabelsvisible=false,
        yticksvisible=false,
        #ylabelfont="Computer Modern Serif"
        #xlabelvisible=false,
        #ylabelvisible=false,
    ),
    Axis3=(
        xgridvisible=false,
        ygridvisible=false,
        zgridvisible=false,
        xspinesvisible=false,
        yspinesvisible=false,
        zspinesvisible=false,
        xticklabelsvisible=false,
        yticklabelsvisible=false,
        zticklabelsvisible=false,
        xticksvisible=false,
        yticksvisible=false,
        zticksvisible=false,
        xlabelvisible=false,
        ylabelvisible=false,
        zlabelvisible=false,
    ),
    #backgroundcolor = GLMk.colorant"rgba(10%, 10%, 10%, 1.0)"
)

GLMk.set_theme!(merge(curr_theme, GLMk.theme_latexfonts()))
#GLMk.set_theme!()#GLMk.themes.dark)

red_cb = GLMk.colorant"rgba(216, 27, 96, 1.0)";
blue_cb = GLMk.colorant"rgba(30, 136, 229, 0.8)";
orange_cb = GLMk.colorant"rgba(206,78,0, 1.0)";
orange2_cb = GLMk.colorant"rgba(255, 193, 7, 1.0)";
green_cb = colorant"rgba(0, 77, 64, 1.0)";

c1 = colorant"rgba(100%, 100%, 100%, 0.0)";
c2 = orange_cb;
c3 = blue_cb;
c4 = red_cb;
oranges = range(c1, stop=c2, length=100);
blues = range(c1, stop=c3, length=100);
reds = range(c1, stop=c4, length=100);

function double_power(r, R, p1, p2)
    return (r / R)^p1 / (1 + (r / R)^(p1 + p2))
end

in_out_colorscheme(outcolor, incolor, eps=0.00000001) =
    GLMk.cgrad([outcolor, incolor, incolor, outcolor], [0.0, eps, 1.0 - eps, 1.0])
out_blue = HSLA(200, 0.9, 0.8, 1.0);
in_blue = HSLA(200, 0.9, 0.4, 1.0);
blue_in_blue = in_out_colorscheme(out_blue, in_blue);

figure = GLMk.Figure(figure_padding = -30, resolution=(500,500*(13.5/10)));
obs_el = π/40
obs_az = -1.08*π
ax = GLMk.Axis3(figure[1, 1], elevation=obs_el, azimuth=π, aspect=(1, 1, 13.5/10))
GLMk.xlims!(ax, (-5, 5))
GLMk.ylims!(ax, (-5, 5))
GLMk.zlims!(ax, (-4.5, 9))

# Blackhole
# Create vertices for a Sphere
r = 1.0+√(1-spin^2)
n = 30
θ = LinRange(0, pi, n)
φ2 = LinRange(0, 2pi, 2 * n)
x2 = [r * cos(φv) * sin(θv) for θv in θ, φv in φ2]
y2 = [r * sin(φv) * sin(θv) for θv in θ, φv in φ2]
z2 = [r * cos(θv) for θv in θ, φv in 2φ2] #.+ 0.2
points = vec([GeometryBasics.Point3f(xv, yv, zv) for (xv, yv, zv) in zip(x2, y2, z2)])

# The coordinates form a matrix, so to connect neighboring vertices with a face
# we can just use the faces of a rectangle with the same dimension as the matrix:
faces = decompose(QuadFace{GLIndex}, Tesselation(Rect(0, 0, 1, 1), size(z2)))

gb_mesh = GeometryBasics.Mesh(points, faces)
GLMk.mesh!(ax, gb_mesh, color=colorant"rgba(10%,10%,10%,1.0)", shading=true, overdraw=false)

# Cone 1
θ = LinRange(0, 0.5, 100)
φ = LinRange(0, 2, 100)
x = sqrt(2) .* 8 .* [cospi(φ) * sinpi(θ) for θ in θ, φ in φ]
y = sqrt(2) .* 8 .* [sinpi(φ) * sinpi(θ) for θ in θ, φ in φ]
z1 = sqrt.(x .^ 2 + y .^ 2) #.+ 0.2
R = 3.0
p1 = 3.0
p2 = 5.0
colorvals1 = double_power.(sqrt.((x .^ 2) .+ (y .^ 2) .+ (z1 .^ 2)), R, p1, p2) 

cone1 = GLMk.surface!(ax, x, y, z1, color=colorvals1, colormap=blues, shading=false, overdraw=false, transparency=true)


# Cone 2
z2 = .-sqrt.(x .^ 2 + y .^ 2)# .+ 0.2
colorvals2 = double_power.(sqrt.((x .^ 2) .+ (y .^ 2) .+ (z2 .^ 2)), R, p1, p2) 
cone2 = GLMk.surface!(ax, x, y, z2, color=colorvals2, colormap=blues, shading=false, overdraw=false, transparency=true)

# Magnetic Field Lines
maglines = load(joinpath((@__DIR__) , "Magnetic_Field_OBJ", "Magnetic Field Highpoly OBJ.obj"))
#maglines = load((@__DIR__) * "/Magnetic_Field_OBJ/maglines1.obj")

rot = Rotations.AngleAxis(π / 2, 1.0, 0.0, 0.0)
#rot = Rotations.AngleAxis(0.0, 1.0, 0.0, 0.0)
rot2 = Rotations.AngleAxis(pi/4, 0.0, 0.0, 1.0)
points = (Ref(rot2) .* (Ref(rot) .* maglines.position)).* 0.2
faces = getfield(getfield(maglines, :simplices), :faces)
mag_mesh = GeometryBasics.Mesh(points, faces)

#GLMk.mesh!(ax, mag_mesh, color=red_cb, linewidth=0.01, shading=false, overdraw=false)
GLMk.mesh!(
    mag_mesh,
    color=[parse(GLMk.Colorant, "rgba(84%, 11%, 38%,$(0.75>(tri[1]/√(tri[1]^2 +tri[2]^2)) > 0.25 ? 1.0 : 0.0))") for tri in mag_mesh.position],
    #color=[parse(GLMk.Colorant, "rgba(84%, 11%, 38%,1)") for tri in mag_mesh.position],
    colormap = :blues, transparency=true,
    shading=false
)
maglines_plot = Ref(ax.scene.plots[end])
#GLMk.rotate!(ax.scene.plots[end],(0,0,1),π/10)



# Intensity profile
#function double_power(r, R, p1, p2)
#    return (r/R)^p1 / (1 + (r/R)^(p1 + p2))
#end
#
#xvals = 0:0.1:12
#zvals = 0:0.1:12
#rvals = sqrt.(xvals.^2 .+ zvals.^2)
#yvals = -10 .* double_power.(rvals, R, p1, p2)
#
#GLMk.lines!(ax, xvals, yvals, zvals, color=:black, linewidth=2)

# Geodesic 1

thickness=0.2
model = Krang.Kerr(spin);
θo = π/10;
rmin = Krang.horizon(model)
#α = sqrt(27/2) + 0.8#1.719
#β = sqrt(27/2) + 0.1
α =-3.874874874874875 #-sqrt(27/2) + 1.9#1.719
β =-1.8008008008008007 #sqrt(27/2) + 0.145

function raytrace(spin, α, β, θo)
    pix = Krang.SlowLightIntensityPixel(Krang.Kerr(spin), α, β, θo)
    return Krang.raytrace.(Ref(pix),range(0.095,1.28,length=50))
end
outvals = raytrace(spin, α, β, θo)
#outvals = Krang.raytrace.(Ref(model), α, β, θo, range(0.095,1.28,length=50))

rvals = []
θvals = []
ϕvals = []
for i in 1:length(outvals)
    push!(rvals, outvals[i][2])
    push!(θvals, outvals[i][3] )
    push!(ϕvals, outvals[i][4])
end

xvals = collect(rvals .* sin.(θvals) .* cos.(ϕvals))
yvals = collect(rvals .* sin.(θvals) .* sin.(ϕvals))
zvals = collect(rvals .* cos.(θvals))

obs_dir = [cos(obs_el)*cos(obs_az), cos(obs_el)*sin(obs_az), sin(obs_el)]
north_dir = LinearAlgebra.normalize([0.,0.,1.0])# .- (LinearAlgebra.dot(obs_dir,[0.,0.,1.0]) .* obs_dir))
dxvals = xvals[begin+1:end] .- xvals[begin:end-1]
dyvals = yvals[begin+1:end] .- yvals[begin:end-1]
dzvals = zvals[begin+1:end] .- zvals[begin:end-1]

pvals = collect.(zip(dxvals, dyvals, dzvals))

euclid_norm(vec) = sum(vec .^ 2)
tick_dir = LinearAlgebra.normalize.(Ref(north_dir) .- (LinearAlgebra.dot.(pvals, Ref(north_dir)) .* pvals ./ euclid_norm.(pvals))) ./ 2# -LinearAlgebra.normalize.((LinearAlgebra.dot.(pvals, Ref(north_dir)) .* pvals)) ./ 2

x_tick_vals = [i[1] for i in tick_dir]
y_tick_vals = [i[2] for i in tick_dir]
z_tick_vals = [i[3] for i in tick_dir]

θvals = acos.([i[3] for i in tick_dir])

s_div=13

#nullhead = GeometryBasics.Mesh([GeometryBasics.Point(0f0,0f0,0f0),], [GeometryBasics.TriangleFace(OffsetInteger{1, UInt32}(0x00000000),OffsetInteger{1, UInt32}(0x00000001),OffsetInteger{1, UInt32}(0x00000002)),])
nullhead = GeometryBasics.Mesh([GeometryBasics.Point(0f0,0f0,0f0),], [GeometryBasics.TriangleFace(OffsetInteger{1, UInt32}(0x00000001),OffsetInteger{1, UInt32}(0x00000002),OffsetInteger{1, UInt32}(0x00000003)),])
GLMk.arrows!(
    ax, 
    xvals[begin+1:2:end-1] .+ (x_tick_vals[begin+1:2:end] ./2), 
    yvals[begin+1:2:end-1] .+ (y_tick_vals[begin+1:2:end] ./2), 
    zvals[begin+1:2:end-1] .+ (z_tick_vals[begin+1:2:end] ./2), 
    x_tick_vals[begin+1:2:end], 
    y_tick_vals[begin+1:2:end], 
    z_tick_vals[begin+1:2:end], 
    color=orange_cb, 
    align=:center, 
    arrowhead=nullhead, 
    shading=true,
)
GLMk.arrows!(
    ax, 
    xvals[begin:end-2], 
    yvals[begin:end-2], 
    zvals[begin:end-2], 
    dxvals[begin:end-1], 
    dyvals[begin:end-1], 
    dzvals[begin:end-1], 
    color=orange_cb, 
    arrowhead=nullhead, 
    shading=true,
    arrowsize=0.25,
)
GLMk.arrows!(
    ax, 
    [xvals[1],], 
    [yvals[1],], 
    [zvals[1],], 
    .-[0.225dxvals[1],], 
    .-[0.225dyvals[1],], 
    .-[0.225dzvals[1],], 
    color=orange_cb, 
    shading=true,
    arrowsize=0.25,
)
GLMk.scatter!(ax, xvals[end], yvals[end], zvals[end], markersize=20, color=orange_cb, depth_shift=-1f0, strokewidth=5, strokecolor=orange2_cb, overdraw=false)
GLMk.scatter!(ax, xvals[begin]-0.35dxvals[1], yvals[begin]-0.35dyvals[1], zvals[begin]-0.35dzvals[1], markersize=0.01, color=orange_cb, depth_shift=-1f0, marker=:x, markerspace=:data, transform_marker=true, rotations=GLMk.Quaternion(3,0,9,0))#, strokewidth=1, strokecolor=:black, overdraw=false)
geodesic1_plots = Ref(ax.scene.plots[end-4:end])


# Geodesic 2
model = Krang.Kerr(spin);
θo = π/10;
rmin = Krang.horizon(model)
ρ = sqrt(27) +0.06
φ = 1.5π/4
α = 4.0990990990990996#cos(φ)*ρ
β = 3.4824824824824825#sin(φ)*ρ

#outvals = Krang.raytrace.(Ref(model), α, β, θo, range(0.085,1.36,length=50))
outvals = raytrace(spin, α, β, θo)


rvals = []
θvals = []
ϕvals = []
for i in 1:length(outvals)
    push!(rvals, outvals[i][2])
    push!(θvals, outvals[i][3] )
    push!(ϕvals, outvals[i][4])
end

xvals = filter(x->!isnan(x),collect(rvals .* sin.(θvals) .* cos.(ϕvals)))
yvals = filter(x->!isnan(x),collect(rvals .* sin.(θvals) .* sin.(ϕvals)))
zvals = filter(x->!isnan(x),collect(rvals .* cos.(θvals)))

obs_dir = [cos(obs_el)*cos(obs_az), cos(obs_el)*sin(obs_az), sin(obs_el)]
north_dir = LinearAlgebra.normalize([0.,0.,1.0] .- (LinearAlgebra.dot(obs_dir,[0.,0.,1.0]) .* obs_dir))
dxvals = xvals[begin+1:end] .- xvals[begin:end-1]
dyvals = yvals[begin+1:end] .- yvals[begin:end-1]
dzvals = zvals[begin+1:end] .- zvals[begin:end-1]

pvals = collect.(zip(dxvals, dyvals, dzvals))

tick_dir = LinearAlgebra.normalize.(Ref(north_dir) .- (LinearAlgebra.dot.(pvals, Ref(north_dir)) .* pvals ./ euclid_norm.(pvals))) ./ 2
x_tick_vals = [i[1] for i in tick_dir]
y_tick_vals = [i[2] for i in tick_dir]
z_tick_vals = [i[3] for i in tick_dir]

θvals = acos.([i[3] for i in tick_dir])

s_div=13

#nullhead = GeometryBasics.Mesh([GeometryBasics.Point(0f0,0f0,0f0),], [GeometryBasics.TriangleFace(OffsetInteger{1, UInt32}(0x00000000),OffsetInteger{1, UInt32}(0x00000001),OffsetInteger{1, UInt32}(0x00000002)),])
GLMk.arrows!(
    ax, 
    xvals[begin+1:2:end-1] .+ (x_tick_vals[begin+1:2:end] ./2), 
    yvals[begin+1:2:end-1] .+ (y_tick_vals[begin+1:2:end] ./2), 
    zvals[begin+1:2:end-1] .+ (z_tick_vals[begin+1:2:end] ./2), 
    x_tick_vals[begin+1:2:end], 
    y_tick_vals[begin+1:2:end], 
    z_tick_vals[begin+1:2:end], 
    color=orange2_cb, 
    align=:center, 
    arrowhead=nullhead, 
    shading=true,
)

GLMk.arrows!(
    ax, 
    xvals[begin:end-2], 
    yvals[begin:end-2], 
    zvals[begin:end-2], 
    dxvals[begin:end-1], 
    dyvals[begin:end-1], 
    dzvals[begin:end-1], 
    color=orange2_cb, 
    arrowhead=nullhead, 
    shading=true,
    arrowsize=0.25,
)
GLMk.arrows!(
    ax, 
    [xvals[1],], 
    [yvals[1],], 
    [zvals[1],], 
    .-[0.25dxvals[1],], 
    .-[0.25dyvals[1],], 
    .-[0.25dzvals[1],], 
    color=orange2_cb, 
    shading=true,
    arrowsize=0.25,
)
GLMk.scatter!(ax, xvals[begin]-0.35dxvals[1], yvals[begin]-0.35dyvals[1], zvals[begin]-0.35dzvals[1], markersize=0.01, color=orange2_cb, depth_shift=-1f0, marker=:x, markerspace=:data, transform_marker=true, rotations=GLMk.Quaternion(3,0,9,0))#, strokewidth=1, strokecolor=:black, overdraw=false)

geodesic2_plots = Ref(ax.scene.plots[end-3:end])

#GLMk.scatter!(ax, xvals[end], yvals[end], zvals[end], markersize=20, color=orange_cb)#, strokewidth=1, strokecolor=:black, overdraw=false)

# Grid Lines
for i in (1.8:0.4:2.6)
    xvals = Vector{Float32}()
    yvals = Vector{Float32}()
    zvals = Vector{Float32}()
    dxvals = Vector{Float32}()
    dyvals = Vector{Float32}()
    dzvals = Vector{Float32}()
    for j in range(0, 2, length=100)
        push!(xvals, i * cospi(j) )
        push!(yvals, i * sinpi(j) )
        push!(zvals, i )#+ 0.2) 
    end

    #GLMk.lines!(
    #    ax, 
    #    xvals, 
    #    yvals, 
    #    zvals, 
    #    color=GLMk.colorant"rgba(10%, 10%, 10%, 0.6)",
    #    linewidth=1
    #)
end

for i in (1.8:0.4:2.6)
    xvals = Vector{Float32}()
    yvals = Vector{Float32}()
    zvals = Vector{Float32}()
    dxvals = Vector{Float32}()
    dyvals = Vector{Float32}()
    dzvals = Vector{Float32}()
    for j in range(0, 2, length=100)
        push!(xvals, i * cospi(j) )
        push!(yvals, i * sinpi(j) )
        push!(zvals, -i)# + 0.2) 
    end

    #GLMk.lines!(
    #    ax, 
    #    xvals, 
    #    yvals, 
    #    zvals, 
    #    color=GLMk.colorant"rgba(10%, 10%, 10%, 0.6)", 
    #    linewidth=2
    #)
end

# Fluid Velocity
xvals = []
yvals = []
zvals = []
dxvals = []
dyvals = []
dzvals = []
for i in 1.8:0.4:2.6
    push!.(Ref(xvals ), i .* cospi.(obs_az/π .+ range(0, 2, length=10)) |> collect)
    push!.(Ref(yvals ), i .* sinpi.(obs_az/π .+ range(0, 2, length=10)) |> collect)
    push!.(Ref(zvals ), i .* ones(10)) 
                     
    push!.(Ref(dxvals), 0.2 .* (cospi.(obs_az/π .+ range(0, 2, length=10) .+ 0.5) .+  cospi.(obs_az/π .+ range(0, 2, length=10))) |> collect)
    push!.(Ref(dyvals), 0.2 .* (sinpi.(obs_az/π .+ range(0, 2, length=10) .+ 0.5) .+  sinpi.(obs_az/π .+ range(0, 2, length=10))) |> collect)
    push!.(Ref(dzvals), 0.2 .* ones(10))

end

for i in 1.8:0.4:2.6
    push!.(Ref(xvals ), .-i .* cospi.(-obs_az/π .+ range(0, 2, length=10)) |> collect)
    push!.(Ref(yvals ), .-i .* sinpi.(-obs_az/π .+ range(0, 2, length=10)) |> collect)
    push!.(Ref(zvals ), .-i .* ones(10))
                     
    push!.(Ref(dxvals), .-0.2 .* (cospi.(-obs_az/π .+ range(0, 2, length=10) .+ 0.5) .+  cospi.(-obs_az/π .+ range(0, 2, length=10))) |> collect)
    push!.(Ref(dyvals), .-0.2 .* (sinpi.(-obs_az/π .+ range(0, 2, length=10) .+ 0.5) .+  sinpi.(-obs_az/π .+ range(0, 2, length=10))) |> collect)
    push!.(Ref(dzvals), .-0.2 .* ones(10) )

end

GLMk.arrows!(
    ax, 
    xvals, 
    yvals, 
    zvals, #.+ 0.2, 
    dxvals, 
    dyvals, 
    dzvals, 
    color=green_cb, 
    shading=false, 
    arrowsize=0.15, 
    align=:tailend
)

# Spin Axis
GLMk.arrows!(
    ax, 
    [0.0,], 
    [0.0,], 
    [0.0,], 
    [0.0,], 
    [0.0,], 
    [4.75,], 
    color=:black, 
    shading=false,
    arrowsize=0.25,
)

GLMk.lines!(
    ax, 
    cospi.(range(0.1, 1.2, length=20)) |> collect, 
    sinpi.(range(0.1, 1.2, length=20)) |> collect, 
    4.5 .* ones(20), #.+ 0.1.*cospi.(range(0.1, 1.2, length=20)) |> collect, 
    color=:black, 
    linewidth=3,
)
GLMk.arrows!(
    ax, 
    [cospi(0.1),], 
    [sinpi(0.1),], 
    [4.5,],# + 0.1.*cospi(0.1),],
    [-0.1*sinpi(0.1),], 
    [-0.1*cos(0.1),], 
    [0.0,],
    color=:black, 
    arrowsize=0.15,
    shading=false
)

GLMk.lines!(
    ax, 
    zeros(10) |> collect,
    4.15 .* cospi.(range(1/2, 3/4, length=10)) |> collect, 
    4.15 .* sinpi.(range(1/2, 3/4, length=10)) |> collect, 
    color=:black, 
    linewidth=2,
)

# Text labels
GLMk.text!(ax, 0.0, 0.2, 5.0;text=GLMk.L"$a$", color=:black, fontsize=30)
GLMk.text!(ax, 1.5r*cos(obs_az)+0.05sin(obs_az), 1.5r*sin(obs_az)-0.05cos(obs_az), -0.40;text=GLMk.L"$M$", color=:white, fontsize=40)
GLMk.text!(ax, 1.3sin(obs_az)+0.5cos(obs_az), 1.3cos(obs_az)+0.5sin(obs_az), 3.2;text=GLMk.L"$\theta_s$", color=:black, fontsize=30, depth_shift=-1f0)
GLMk.text!(ax, 4.5sin(obs_az), 4.5cos(obs_az), -0.5;text=GLMk.L"$\vec{B}$", color=:black, fontsize=30)

GLMk.lines!(
    ax, 
    -2ones(20),
    2.4 .+ range(0, 0.5, length=20),
    5.8 .+ 0.1 .* sinpi.(range(0.5, 2.5, length=20)) |> collect, 
    color=:black, 
    linewidth=2,
)

GLMk.arrows!(
    ax, 
    [-2,],
    [2.4, ],
    [5.8 + 0.1sinpi(0.5),],
    [0.0,],
    [-0.01,],
    [0.0],
    color=:black, 
    linewidth=0.1,
)
GLMk.text!(ax, -2, 3.5, 5.6;text=GLMk.L"$\vec{P}$", color=:black, fontsize=30)

GLMk.lines!(
    ax, 
    -3.0 .* ones(20),
    0.1 .* sinpi.(range(0.5, 2.5, length=20)) |> collect, 
    .-3.5 .+ range(0, 0.5, length=20),
    color=:black, 
    linewidth=2,
)

GLMk.arrows!(
    ax, 
    [-3.0,],
    [0.1 .* sinpi.(0.5), ],
    [-3.0,],
    [0.0,],
    [0.0,],
    [0.01],
    color=:black, 
    linewidth=0.1,
)
GLMk.text!(
    ax, 
    -3,
    0.2 .+ 0.1 .* sinpi.(0.5),
    -4.5;
    text=GLMk.L"$\vec{k}$", 
    color=:black, 
    fontsize=30
)

# Emission profile envelope
r = range(0, 10, length=200)

#ax = GLMk.Axis(
#    figure[2, 1], 
#    aspect=5, 
#    xlabel=GLMk.L"$r$", 
#    ylabel=GLMk.L"$\mathcal{J}$",
#    xticks=([R,], [GLMk.L"R",]),
#    title=GLMk.L"\text{Enveloping Intensity Profile}",
#    topspinevisible=false,
#    rightspinevisible=false,
#    )
#prof = double_power.(r, R, p1, p2) 
#GLMk.lines!(ax, r, prof, color=blue_cb, linewidth=2, xscale=[0.5])
#
#GLMk.text!(ax, 5, 0.05;text=GLMk.L"\mathcal{J}= \frac{(r/R)^{p_1}}{1+(r/R)^{p_1+p_2}}", fontsize=20)
#GLMk.rowgap!(figure.layout, 1, GLMk.Fixed(0))
#GLMk.rowsize!(figure.layout, 1, GLMk.Auto(5))



using Pyehtim
img = ehtim.image.load_fits(joinpath((@__DIR__), "best.fits")).regrid_image(80*ehtim.RADPERUAS, 1000)
#img.display()
intensity = reshape(pyconvert(Vector{Float64},img.imvec), (1000,1000))
#ax = GLMk.Axis3(
#    figure[1,1], 
#    #azimuth=0π/180,
#    #elevation=4.25π/4,
#    #aspect=1, 
#    #width=GLMk.Relative(0.1), 
#    #height=GLMk.Relative(0.1), 
#    #halign=1.50, 
#    #valign=-4.00, 
#    #backgroundcolor=:white,
#    )
#GLMk.empty!(figure)
#GLMk.hidedecorations!(ax)
#ax = GLMk.LScene(figure[1,1], show_axis=false)
curr_scene = ax.scene
#cam = GLMk.Camera3D(curr_scene)#, projectiontype=0)
#
function mat_to_quat(m)
    if (m[3,3] < 0) 
        if (m[1,1] >m[2,2]) 
            t = 1 + m[1,1] -m[2,2] -m[3,3];
            q = GLMk.Quaternion( t, m[1,2]+m[2,1], m[3,1]+m[1,3], m[2,3]-m[3,2] );
        
        else 
            t = 1 -m[1,1] + m[2,2] -m[3,3];
            q = GLMk.Quaternion( m[1,2]+m[2,1], t, m[2,3]+m[3,2], m[3,1]-m[1,3] );
        end
    
    else 
        if (m[1,1] < -m[2,2]) 
            t = 1 -m[1,1] -m[2,2] + m[3,3];
            q = GLMk.Quaternion( m[3,1]+m[1,3], m[2,3]+m[3,2], t, m[1,2]-m[2,1] );
        
        else 
            t = 1 + m[1,1] + m[2,2] + m[3,3];
            q = GLMk.Quaternion( m[2,3]-m[3,2], m[3,1]-m[1,3], m[1,2]-m[2,1], t );
        end
    end
    q = GLMk.Quaternion((q.data .* (0.5 / sqrt(t))...))
end
m(x,y,z)=
[
    cos(z)*cos(y) -cos(y)*sin(z) -sin(y)
    cos(z)*cos(x)*sin(y)+sin(z)*sin(x) -cos(x)*sin(z)*sin(y)+cos(z)*sin(x) cos(y)*cos(x) 
    -cos(x)*sin(z)+cos(z)*sin(y)*sin(x) -cos(z)*cos(x)-sin(z)*sin(y)*sin(x) cos(y)*sin(x)
]
m(z,y) = [
    cos(z)*cos(y) -cos(y)*sin(z) -sin(y) 
    cos(z)*sin(y) -sin(z)*sin(y) cos(y) 
    -sin(z) -cos(z) 0
]
mat_to_quat(m(0,0))

hm = GLMk.heatmap!(ax,range(-5,5,length=400),range(-5,5,length=400), intensity; colormap=:afmhot)#, show_axis = false)
#GLMk.arrows!(ax, [4.0,4.0],[3.25,3.25], zeros(2), [0,-2], [-2, 0], zeros(2), color=:white, linewidth=0.1, arrowsize=0.1, depth_shift=-1f0, shading=false)#, arrowhead=nullhead)
#GLMk.text!(ax, 1.5, 3.0, 0; text=GLMk.L"$\beta$", color=:white, fontsize=0.2, depth_shift=-1f0, markerspace=:data, rotation=mat_to_quat(m(3π/2-π/10,π,π/2)))
#GLMk.text!(ax, 4.5, 0.5, 0; text=GLMk.L"$\alpha$", color=:white, fontsize=0.2, markerspace=:data, depth_shift=-1f0, rotation=mat_to_quat(m(3π/2-π/10,π,π/2)))

#cam.lookat[]=[0,0,0]
#cam.eyeposition[] = [0,1,10]
#GLMk.update_cam!(ax.scene)
GLMk.scale!.(curr_scene.plots[end:end],1.4, 1.4, 1.4)
#GLMk.rotate!.(curr_scene.plots[end-3:end],Ref((0,1,0)),Ref(π/2))
GLMk.rotate!.(curr_scene.plots[end:end],Ref((0,1,0)),Ref(π/10))
GLMk.translate!.(curr_scene.plots[end:end],Ref((0.75,0,8)))
#GLMk.rotate!.(curr_scene.plots,Ref((0,0,1)),Ref(-π/10))
#GLMk.rotate!(maglines_plot.x,(0,0,1),-π/4)
GLMk.rotate!.(geodesic1_plots.x,Ref((0,0,1)),Ref(2π/10))
GLMk.rotate!.(geodesic2_plots.x,Ref((0,0,1)),Ref(2π/10))
GLMk.scale!.(geodesic1_plots.x,0.7,0.7,0.7)
GLMk.scale!.(geodesic2_plots.x,0.7,0.7,0.7)

display(figure)



save(joinpath((@__DIR__), "pedagogical_figure2.png"), figure, pt_per_unit = 10.)

