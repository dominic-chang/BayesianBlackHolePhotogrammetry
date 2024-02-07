export GPUModel
using Metal
using VLBISkyModels

"""
    GPUModel
Experimental model wrapper than enables multi-threading when evaluating `intensitymap`
"""
struct GPUModel{M} <: ComradeBase.AbstractModel
    model::M
end

@inline GPUModel(model::GPUModel) = model
@inline GPUModel(model::M) where {M<:VLBISkyModels.CompositeModel} = GPUModel(ComradeBase.imanalytic(M), model)
@inline GPUModel(::ComradeBase.IsAnalytic, model::VLBISkyModels.CompositeModel)        = GPUModel{typeof(model)}(model)
@inline GPUModel(::ComradeBase.NotAnalytic, model::VLBISkyModels.AddModel)             = AddModel(GPUModel(model.m1), GPUModel(model.m2))
@inline GPUModel(::ComradeBase.NotAnalytic, model::VLBISkyModels.ConvolvedModel)       = ConvolvedModel(GPUModel(model.m1), GPUModel(model.m2))
@inline GPUModel(model::VLBISkyModels.ModelImage) = VLBISkyModels.@set model.model = GPUModel(model.model)

Base.@constprop :aggressive @inline visanalytic(::Type{<:GPUModel{M}}) where {M} = ComradeBase.visanalytic(M)
Base.@constprop :aggressive @inline imanalytic(::Type{<:GPUModel{M}}) where {M} = ComradeBase.imanalytic(M)
Base.@constprop :aggressive @inline ispolarized(::Type{<:GPUModel{M}}) where {M} = ComradeBase.ispolarized(M)

@inline function ComradeBase.visibility_point(m::GPUModel, u, v, time, freq) ComradeBase.visibility_point(m.model, u, v, time, freq) end
@inline function ComradeBase.intensity_point(m::GPUModel, p)  ComradeBase.intensity_point(m.model, p) end

@inline ComradeBase.radialextent(m::GPUModel) = ComradeBase.radialextent(basemodel(m))
@inline ComradeBase.flux(m::GPUModel) = ComradeBase.flux(ComradeBase.basemodel(m))

function ComradeBase.intensitymap_analytic(s::GPUModel, g::RectiGrid)
    T = typeof(ComradeBase.intensity_point(s, (X=g.X[begin], Y=g.Y[begin])))
    img = ComradeBase.IntensityMap(Array{T}(undef, length(g.X), length(g.Y)), g)
    return ComradeBase.intensitymap_analytic!(img, s)
end

ComradeBase.intensitymap(m, p, threaded::Bool) = intensitymap(m, p, static(threaded))
ComradeBase.intensitymap(m, p, ::VLBISkyModels.False) = intensitymap(m, p)
ComradeBase.intensitymap(m, p, ::VLBISkyModels.True)  = intensitymap(GPUModel(m), p)

ComradeBase.intensitymap!(img::IntensityMapTypes, m, threaded::Bool) = intensitymap!(img, m, static(threaded))
ComradeBase.intensitymap!(img::IntensityMapTypes, m, ::VLBISkyModels.False) = intensitymap!(img, m)
ComradeBase.intensitymap!(img::IntensityMapTypes, m, ::VLBISkyModels.True)  = intensitymap!(img, GPUModel(m))

function ComradeBase.intensitymap!(::ComradeBase.IsAnalytic, img::ComradeBase.IntensityMap, s::GPUModel)
    dx, dy = pixelsizes(img)
    mm = Base.Fix1(ComradeBase.intensity_point, s)
    g = MtlArray(imagegrid(img))
    img .= mm.(g)*dx*dy
    println("here")
    return img
end

function fouriermap(::ComradeBase.IsAnalytic, m::GPUModel, dims::ComradeBase.AbstractGrid)
    X = dims.X
    Y = dims.Y
    uu,vv = uviterator(length(X), step(X), length(Y), step(Y))
    # uvgrid = ComradeBase.grid(U=uu, V=vv)
    T = typeof(visibility_point(m, uu[1], vv[1], 0, 0))
    vis = similar(uu, T, (length(uu), length(vv)))
    Threads.@threads for I in CartesianIndices(vis)
        ix, iy = Tuple(I)
        vis[I] = visibility_point(m, uu[ix], vv[iy], 0, 0)
    end
    return IntensityMap(vis, (U=uu, V=vv))
end