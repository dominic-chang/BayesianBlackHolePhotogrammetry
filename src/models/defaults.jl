abstract type AbstractModelInfo end
struct ModelInfo{M,T,F} <: AbstractModelInfo
    model::M
    name::String
    lower::NamedTuple{T}
    upper::NamedTuple{T}
    file::String
    blur::F
    function ModelInfo(model::M, name::String, lower::NamedTuple{T}, upper::NamedTuple{T}, file::String, blur::F=zero(F)) where {M,T,F}
        check = match(r"i[0-9]{2}", file)

        if !isnothing(check)
            inc = parse(Int, check.match[2:end])
            lower = NamedTuple{keys(lower)}([i == :θo ? max((inc - 20.0f0) / 180.0f0 * π, 0.0f0) : lower[i] for i in keys(lower)])
            upper = NamedTuple{keys(upper)}([i == :θo ? min((inc + 20.0f0) / 180.0f0 * π, π/2f0) : upper[i] for i in keys(upper)])
        end
        new{M,T,F}(model, name, lower, upper, file, blur)
    end
end

lower_eq = (
    m_d=1.5f0,
    spin=-1.0f0,
    θo=0.0f0,
    rpeak=1f0,
    p1=0.1f0,
    p2=1.0f0,
    χ=-1.0f0π,
    ι=Float32(0f0),
    βv=0.01f0,
    σ=-1.0f0,
    x0=-2.0f0,
    y0=-2.0f0
)

upper_eq = (
    m_d=(8.0f0),
    spin=-0.01f0,
    θo=40.0f0 / 180 * π,
    rpeak=10f0,
    p1=10.0f0,
    p2=10.0f0,
    χ=1.0f0π,
    ι=π / 2.0f0,
    βv=0.9f0,
    σ=3.0f0,
    x0=2.0f0,
    y0=2.0f0
)

lower_dc = (
    m_d=1.5f0,
    spin=-1.0f0,
    θo=0.0f0,
    θs=20f0/180*π,
    rpeak=1f0,
    p1=0.1f0,
    p2=1.0f0,
    χ=-1.0f0π,
    ι=Float32(0f0),
    βv=0.01f0,
    σ=-1.0f0,
    η=-1f0π,
    x0=-2.0f0,
    y0=-2.0f0
)


upper_dc = (
    m_d=(8.0f0),
    spin=-0.01f0,
    θo=40.0f0 / 180 * π,
    θs=π / 2f0,
    rpeak=10f0,
    p1=10.0f0,
    p2=10.0f0,
    χ=1.0f0π,
    ι=π / 2.0f0,
    βv=0.9f0,
    σ=3.0f0,
    η=1f0π,
    x0=2.0f0,
    y0=2.0f0
)

lower_eq_dc = (
    m_d=1.5f0,
    spin=-1.0f0,
    θo=0.0f0,
    θs=20f0/180*π,
    r_cone=1.0f0,
    p1_cone=0.1f0,
    p2_cone=1.0f0,
    χ_cone=-1.0f0π,
    ι_cone=Float32(0f0),
    βv_cone=0.01f0,
    spec_cone=-1.0f0,
    η_cone=-1f0π,
    r_disk=1.0f0,
    p1_disk=0.1f0,
    p2_disk=1.0f0,
    χ_disk=-1.0f0π,
    ι_disk=Float32(0f0),
    βv_disk=0.01f0,
    spec_disk=-1.0f0,
    η_disk=-1f0π,
    rJ=0.0f0,
    xtrans=-2.0f0,
    ytrans=-2.0f0
)

upper_eq_dc = (
    m_d=8.0f0,
    spin=-0.01f0,
    θo=40.0f0 / 180 * π,
    θs=π / 2f0,
    r_cone=10.0f0,
    p1_cone=10.0f0,
    p2_cone=10.0f0,
    χ_cone=1.0f0π,
    ι_cone=π / 2f0,
    βv_cone=0.9f0,
    spec_cone=3.0f0,
    η_cone=1f0π,
    r_disk=10.0f0,
    p1_disk=10.0f0,
    p2_disk=10.0f0,
    χ_disk=1.0f0π,
    ι_disk=π / 2f0,
    βv_disk=0.9f0,
    spec_disk=3.0f0,
    η_disk=1.0f0π,
    rJ=1.0f0,
    xtrans=2.0f0,
    ytrans=2.0f0
)
