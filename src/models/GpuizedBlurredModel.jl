export GPUBlurModel
using Metal
using VLBISkyModels
using ImageFiltering

"""
    GPUBlurModel
Experimental model wrapper than enables multi-threading when evaluating `intensitymap`
"""
struct GPUBlurModel{M,T} <: ComradeBase.AbstractModel
    model::M
    blur::T
end

@inline GPUBlurModel(model::GPUBlurModel) = model
@inline GPUBlurModel(model::M) where {M<:VLBISkyModels.CompositeModel} = GPUBlurModel(ComradeBase.imanalytic(M), model)
@inline GPUBlurModel(::ComradeBase.IsAnalytic, model::VLBISkyModels.CompositeModel) = GPUBlurModel{typeof(model)}(model)
@inline GPUBlurModel(::ComradeBase.NotAnalytic, model::VLBISkyModels.AddModel) = AddModel(GPUBlurModel(model.m1), GPUBlurModel(model.m2))
@inline GPUBlurModel(::ComradeBase.NotAnalytic, model::VLBISkyModels.ConvolvedModel) = VLBISkyModels.ConvolvedModel(GPUBlurModel(model.m1), GPUBlurModel(model.m2))
@inline GPUBlurModel(model::VLBISkyModels.ModelImage) = VLBISkyModels.@set model.model = GPUBlurModel(model.model)

Base.@constprop :aggressive @inline visanalytic(::Type{<:GPUBlurModel{M}}) where {M} = ComradeBase.visanalytic(M)
Base.@constprop :aggressive @inline imanalytic(::Type{<:GPUBlurModel{M}}) where {M} = ComradeBase.imanalytic(M)
Base.@constprop :aggressive @inline ispolarized(::Type{<:GPUBlurModel{M}}) where {M} = ComradeBase.ispolarized(M)

@inline function ComradeBase.visibility_point(m::GPUBlurModel, u, v, time, freq)
    ComradeBase.visibility_point(m.model, u, v, time, freq)
end
@inline function ComradeBase.intensity_point(m::GPUBlurModel, p)
    ComradeBase.intensity_point(m.model, p)
end

@inline ComradeBase.radialextent(m::GPUBlurModel) = ComradeBase.radialextent(basemodel(m))
@inline ComradeBase.flux(m::GPUBlurModel) = ComradeBase.flux(ComradeBase.basemodel(m))

function ComradeBase.intensitymap_analytic(s::GPUBlurModel, g::RectiGrid)
    T = typeof(ComradeBase.intensity_point(s, (X=g.X[begin], Y=g.Y[begin])))
    img = ComradeBase.IntensityMap(Array{T}(undef, length(g.X), length(g.Y)), g)
    return ComradeBase.intensitymap_analytic!(img, s)
end


struct GPUizeAndBlur end
#ComradeBase.intensitymap(m, p, threaded::Bool) = intensitymap(m, p, static(threaded))
ComradeBase.intensitymap(m::GPUBlurModel, p, threaded) = intensitymap(m, p, GPUized())
#ComradeBase.intensitymap(m, p, ::VLBISkyModels.False) = intensitymap(m, p)
#ComradeBase.intensitymap(m, p, ::VLBISkyModels.True) = intensitymap(GPUBlurModel(m), p)

#ComradeBase.intensitymap!(img::IntensityMapTypes, m, threaded::Bool) = intensitymap!(img, m, static(threaded))
#ComradeBase.intensitymap!(img::IntensityMapTypes, m, ::VLBISkyModels.False) = intensitymap!(img, m)
#ComradeBase.intensitymap!(img::IntensityMapTypes, m, ::VLBISkyModels.True) = intensitymap!(img, ThreadedModel(m))
ComradeBase.intensitymap!(img::IntensityMapTypes, m, ::GPUizeAndBlur) = ComradeBase.intensitymap!(img, GPUBlurModel(m))
ComradeBase.intensitymap!(img::IntensityMapTypes, m::GPUBlurModel) = ComradeBase.intensitymap_analytic!(img, m)

function ComradeBase.intensitymap_analytic!(inimg::ComradeBase.IntensityMap, s::GPUBlurModel)
    dx, dy = pixelsizes(inimg)
    mm = Base.Fix1(ComradeBase.intensity_point, s)
    g = MtlArray(imagegrid(inimg) |> collect)
    inimg .= Array(mm.(g) * dx * dy)
    blur = s.blur
    σ_px = blur / (2 * sqrt(2 * log(2f0))) / abs(dx)
    σ_py = blur / (2 * sqrt(2 * log(2f0))) / abs(dy)

    # Now I need to pick my kernel size. I am going out to 5σ for the
    # gaussian kernel. I have to add one for the convolution to play nice
    nkern = Int(floor(σ_px) * 10 + 1)
    inimg .= imfilter(parent(inimg),
        Kernel.gaussian((σ_py, σ_px), (nkern, nkern)),
        Fill(0.0, inimg),
        Algorithm.FFT()
    )
    return inimg

end

function fouriermap(::ComradeBase.IsAnalytic, m::GPUBlurModel, dims::ComradeBase.AbstractGrid)
    X = dims.X
    Y = dims.Y
    uu, vv = uviterator(length(X), step(X), length(Y), step(Y))
    # uvgrid = ComradeBase.grid(U=uu, V=vv)
    T = typeof(visibility_point(m, uu[1], vv[1], 0, 0))
    vis = similar(uu, T, (length(uu), length(vv)))
    Threads.@threads for I in CartesianIndices(vis)
        ix, iy = Tuple(I)
        vis[I] = visibility_point(m, uu[ix], vv[iy], 0, 0)
    end
    return IntensityMap(vis, (U=uu, V=vv))
end
