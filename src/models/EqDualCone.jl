using Krang
using VLBIImagePriors
using ComradeBase

struct EqDualCone{T, G1, G2, M} <: ComradeBase.AbstractModel
    spin::T
    θo::T
    θs::T
    r_cone::T
    p1_cone::T
    p2_cone::T
    χ_cone::T
    ι_cone::T
    βv_cone::T
    spec_cone::T
    η_cone::T
    r_disk::T
    p1_disk::T
    p2_disk::T
    χ_disk::T
    ι_disk::T
    βv_disk::T
    spec_disk::T
    η_disk::T
    rJ::T
    mesh_cone::Krang.Mesh{G1, M}
    mesh_disk::Krang.Mesh{G2, M}

    function EqDualCone(
        spin::T,
        θo::T,
        θs::T,
        r_cone::T,
        p1_cone::T,
        p2_cone::T,
        χ_cone::T,
        ι_cone::T,
        βv_cone::T,
        spec_cone::T,
        η_cone::T,
        r_disk::T,
        p1_disk::T,
        p2_disk::T,
        χ_disk::T,
        ι_disk::T,
        βv_disk::T,
        spec_disk::T,
        η_disk::T,
        rJ::T,
        nmax::Int
    ) where T
        η2 = π-η_cone
        magfield0 = Krang.SVector(sin(ι_disk)*cos(η_disk), sin(ι_disk)*sin(η_disk), cos(ι_disk));
        magfield1 = Krang.SVector(sin(ι_cone)*cos(η_cone), sin(ι_cone)*sin(η_cone), cos(ι_cone));
        magfield2 = Krang.SVector(sin(ι_cone)*cos(η2), sin(ι_cone)*sin(η2), cos(ι_cone));

        vel_disk = Krang.SVector(βv_disk, T(π/2), χ_disk);
        vel_cone = Krang.SVector(βv_cone, T(π/2), χ_cone);

        material = Krang.ElectronSynchrotronPowerLawIntensity();

        # Create the Geometry
        profile_cone(r) = let  R = r_cone, p1 = p1_cone, p2 = p2_cone
            ((r/R)^p1)/(1+(r/R)^(p1+p2))
        end
        profile_disk(r) = let  R = r_disk, p1 = p1_disk, p2 = p2_disk
            ((r/R)^p1)/(1+(r/R)^(p1+p2))
        end

        subimgs = (i for i in 0:nmax)
        geometry0 = Krang.ConeGeometry(T(π/2), (magfield0, vel_disk, subimgs, profile_disk, spec_disk))
        geometry1 = Krang.ConeGeometry(θs, (magfield1, vel_cone, subimgs, profile_cone, spec_cone))
        geometry2 = Krang.ConeGeometry(π-θs, (magfield2, vel_cone, subimgs, profile_cone, spec_cone))
        geometry = geometry1 ⊕ geometry2

        mesh_cone = Krang.Mesh(geometry, material)
        mesh_disk = Krang.Mesh(geometry0, material)

        return new{T, typeof(geometry), typeof(geometry0), typeof(material)}(
            spin,
            θo,
            θs,
            r_cone,
            p1_cone,
            p2_cone,
            χ_cone,
            ι_cone,
            βv_cone,
            spec_cone,
            η_cone,
            r_disk,
            p1_disk,
            p2_disk,
            χ_disk,
            ι_disk,
            βv_disk,
            spec_disk,
            η_disk,
            rJ,
            mesh_cone,
            mesh_disk
        )
    end
end
function EqDualCone(nmax::Int, θ::NamedTuple)
    (;
        spin,
        θo,
        θs,
        r_cone,
        p1_cone,
        p2_cone,
        χ_cone,
        ι_cone,
        βv_cone,
        spec_cone,
        η_cone,
        r_disk,
        p1_disk,
        p2_disk,
        χ_disk,
        ι_disk,
        βv_disk,
        spec_disk,
        η_disk,
        rJ,
    ) = θ
    return EqDualCone(
        spin,
        θo,
        θs,
        r_cone,
        p1_cone,
        p2_cone,
        χ_cone,
        ι_cone,
        βv_cone,
        spec_cone,
        η_cone,
        r_disk,
        p1_disk,
        p2_disk,
        χ_disk,
        ι_disk,
        βv_disk,
        spec_disk,
        η_disk,
        rJ,
        nmax
    )
end

ComradeBase.visanalytic(::Type{<:EqDualCone}) = ComradeBase.NotAnalytic()
ComradeBase.imanalytic(::Type{<:EqDualCone}) = ComradeBase.IsAnalytic()
ComradeBase.isprimitive(::Type{<:EqDualCone}) = ComradeBase.IsPrimitive()

@inline function ComradeBase.intensity_point(m::EqDualCone{T,G1,G2,M}, p) where {T,G1,G2,M}
    (;X, Y) = p
    (;rJ, mesh_cone, mesh_disk) = m
    pix = Krang.IntensityPixel(Krang.Kerr(m.spin), -X, Y, m.θo)
    return (one(T)-rJ)*mesh_cone.material(pix, mesh_cone.geometry) + rJ*mesh_disk.material(pix, mesh_disk.geometry)

end
