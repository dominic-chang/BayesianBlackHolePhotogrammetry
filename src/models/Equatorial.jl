struct Equatorial{T,G,M} <: ComradeBase.AbstractModel
    spin::T
    θo::T
    rpeak::T
    p1::T
    p2::T
    χ::T
    ι::T
    βv::T
    spec::T
    mesh::Krang.Mesh{G,M}

    function Equatorial(spin::T, θo::T, rpeak::T, p1::T, p2::T, χ::T, ι::T, βv::T, spec::T, nmax::Int) where {T}
        η = π + χ
        magfield = Krang.SVector(sin(ι) * cos(η), sin(ι) * sin(η), cos(ι))
        vel = Krang.SVector(βv, T(π / 2), χ)
        material = Krang.ElectronSynchrotronPowerLawIntensity()

        # Create the Geometry
        profile(r) =
            let R = rpeak, p1 = p1, p2 = p2
                ((r / R)^p1) / (1 + (r / R)^(p1 + p2))
            end
        subimgs = (i for i in 0:nmax)
        geometry = Krang.ConeGeometry(T(π / 2), (magfield, vel, subimgs, profile, spec))
        mesh = Krang.Mesh(geometry, material)

        return new{T,typeof(geometry),typeof(material)}(
            spin,
            θo,
            rpeak,
            p1,
            p2,
            χ,
            ι,
            βv,
            spec,
            mesh
        )
    end
end
function Equatorial(nmax::Int, θ::NamedTuple)
    (;
        spin,
        θo,
        rpeak,
        p1,
        p2,
        χ,
        ι,
        βv,
        spec,
    ) = θ
    return Equatorial(
        spin,
        θo,
        rpeak,
        p1,
        p2,
        χ,
        ι,
        βv,
        spec,
        nmax
    )
end

ComradeBase.visanalytic(::Type{<:Equatorial}) = ComradeBase.NotAnalytic()
ComradeBase.imanalytic(::Type{<:Equatorial}) = ComradeBase.IsAnalytic()
ComradeBase.isprimitive(::Type{<:Equatorial}) = ComradeBase.IsPrimitive()

@inline function ComradeBase.intensity_point(m::Equatorial{T,G,M}, p) where {T,G,M}
    (; X, Y) = p
    (; mesh) = m
    pix = Krang.IntensityPixel(Krang.Kerr(m.spin), -X, Y, m.θo)
    return mesh.material(pix, (mesh.geometry))
end
