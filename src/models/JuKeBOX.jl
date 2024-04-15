struct JuKeBOX{T} <: ComradeBase.AbstractModel
    spin::T
    θo::T
    θs::T
    rpeak::T
    p1::T
    p2::T
    χ::T
    ι::T
    βv::T
    spec::T
    η::T
    nmax::Int
end

function JuKeBOX(nmax::Int, θ::NamedTuple)
    (;
        spin,
        θo,
        θs,
        rpeak,
        p1,
        p2,
        χ,
        ι,
        βv,
        spec,
        η,
    ) = θ
    return JuKeBOX(
        spin,
        θo,
        θs,
        rpeak,
        p1,
        p2,
        χ,
        ι,
        βv,
        spec,
        η,
        nmax
    )
end

ComradeBase.visanalytic(::Type{<:JuKeBOX}) = ComradeBase.NotAnalytic()
ComradeBase.imanalytic(::Type{<:JuKeBOX}) = ComradeBase.IsAnalytic()
ComradeBase.isprimitive(::Type{<:JuKeBOX}) = ComradeBase.IsPrimitive()

@inline function ComradeBase.intensity_point(m::JuKeBOX{T}, p) where {T}
    (; X, Y) = p
    (;ι, η, χ, βv, θs, rpeak, p1, p2, spec, nmax) = m 

    #η2 = π - η
    magfield1 = Krang.SVector(sin(ι) * cos(η), sin(ι) * sin(η), cos(ι))
    magfield2 = Krang.SVector(-sin(ι) * cos(η), -sin(ι) * sin(η), cos(ι))
    vel = Krang.SVector(βv, T(π / 2), χ)
    material = Krang.ElectronSynchrotronPowerLawIntensity()

    # Create the Geometry
    @inline profile(r) = let R=rpeak, p1=p1, p2=p2
        return (r/R)^p1/(1+(r/R)^(p1+p2))
    end
        
    subimgs = (i for i in 0:nmax)
    geometry1 = Krang.ConeGeometry((θs), (magfield1, vel, subimgs, profile, spec))
    geometry2 = Krang.ConeGeometry((π - θs), (magfield2, vel, subimgs, profile, spec))
    geometry = geometry1 ⊕ geometry2

    mesh = Krang.Mesh(geometry, material)

    pix = Krang.IntensityPixel(Krang.Kerr(m.spin), -X, Y, m.θo)
    ans = mesh.material(pix, (mesh.geometry))
    return isnan(ans) ? zero(T) : ans
end

function VLBISkyModels.__extract_tangent(dm::JuKeBOX)
    ntm = VLBISkyModels.NamedTupleTools.ntfromstruct(dm)
    if ntm isa NamedTuple{(), Tuple{}}
        tbm = VLBISkyModels.ChainRulesCore.ZeroTangent()
    else
        tbm = VLBISkyModels.ChainRulesCore.Tangent{typeof(dm)}(;ntm...)
    end
    return tbm
end