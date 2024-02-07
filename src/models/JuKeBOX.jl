using Krang
using VLBIImagePriors

struct JuKeBOX{T, G, M} <: ComradeBase.AbstractModel
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
    mesh::Krang.Mesh{G, M}

    function JuKeBOX(spin::T,θo::T,θs::T,rpeak::T,p1::T,p2::T,χ::T,ι::T,βv::T,spec::T,η::T, nmax::Int) where T
        η2 = π-η
        magfield1 = Krang.SVector(sin(ι)*cos(η), sin(ι)*sin(η), cos(ι));
        magfield2 = Krang.SVector(sin(ι)*cos(η2), sin(ι)*sin(η2), cos(ι));
        vel = Krang.SVector(βv, T(π/2), χ);
        material = Krang.ElectronSynchrotronPowerLawIntensity();

        # Create the Geometry
        profile(r) = let  R = rpeak, p1 = p1, p2 = p2
            ((r/R)^p1)/(1+(r/R)^(p1+p2))
        end
        subimgs = (i for i in 0:nmax)
        geometry1 = Krang.ConeGeometry((θs), (magfield1, vel, subimgs, profile, spec))
        geometry2 = Krang.ConeGeometry((π-θs), (magfield2, vel, subimgs, profile, spec))
        geometry = geometry1 ⊕ geometry2

        mesh = Krang.Mesh(geometry, material)

        return new{T, typeof(geometry), typeof(material)}(
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
            mesh
        )
    end
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

function ComradeBase.intensity_point(m::JuKeBOX{T,G,M}, p) where {T,G,M}
    (;X, Y) = p
    (;f,mesh) = m
    pix = Krang.IntensityPixel(Krang.Kerr(m.spin), -X, Y, m.θo)
    return f*mesh.material(pix, (mesh.geometry))
end
