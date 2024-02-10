using VLBISkyModels
using BlackBoxOptim


function JuKeBOX(m_d::T,spin::T,θo::T,θs::T,rpeak::T,p1::T,p2::T,χ::T,ι::T,βv::T,spec::T,η::T, x0::T,y0::T) where T
    return modify(JuKeBOX(spin,θo,θs,rpeak,p1,p2,χ,ι,βv,spec, η, 2),
    Stretch(m_d, m_d), Shift(x0, y0))
end
function EqDualCone(
        m_d::T,
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
        x0::T,
        y0::T
    ) where T
    return modify(EqDualCone(
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
        2
    ),
    Stretch(m_d, m_d), Shift(x0, y0))
end

function dual_cone_temp(θ)
    return GPUModel(
        JuKeBOX(
            Float32(θ.m_d),
            Float32(θ.spin),
            Float32(θ.θo),
            Float32(θ.θs),
            Float32(θ.rpeak),
            Float32(θ.p1),
            Float32(θ.p2),
            Float32(θ.χ),
            Float32(θ.ι),
            Float32(θ.βv),
            Float32(θ.σ),
            Float32(θ.η),
            Float32(θ.x0),
            Float32(θ.y0)
        )
    )
end

function eq_dual_cone_temp(θ)
    return GPUModel(
        EqDualCone(
            Float32(θ.m_d),
            Float32(θ.spin),
            Float32(θ.θo),
            Float32(θ.θs),
            Float32(θ.r_cone),
            Float32(θ.p1_cone),
            Float32(θ.p2_cone),
            Float32(θ.χ_cone),
            Float32(θ.ι_cone),
            Float32(θ.βv_cone),
            Float32(θ.spec_cone),
            Float32(θ.η_cone),
            Float32(θ.r_disk),
            Float32(θ.p1_disk),
            Float32(θ.p2_disk),
            Float32(θ.χ_disk),
            Float32(θ.ι_disk),
            Float32(θ.βv_disk),
            Float32(θ.spec_disk),
            Float32(θ.η_disk),
            Float32(θ.rJ),
            Float32(θ.xtrans),
            Float32(θ.ytrans)
        )
    )
end

function create_dual_cone_blur_model(blur)
    return θ -> GPUBlurModel(
        JuKeBOX(
            Float32(θ.m_d),
            Float32(θ.spin),
            Float32(θ.θo),
            Float32(θ.θs),
            Float32(θ.rpeak),
            Float32(θ.p1),
            Float32(θ.p2),
            Float32(θ.χ),
            Float32(θ.ι),
            Float32(θ.βv),
            Float32(θ.σ),
            Float32(θ.η),
            Float32(θ.x0),
            Float32(θ.y0)
        ), blur
    )
end

function create_eq_dual_cone_blur_model(blur)
    return θ -> GPUBlurModel(
        EqDualCone(
            Float32(θ.m_d),
            Float32(θ.spin),
            Float32(θ.θo),
            Float32(θ.θs),
            Float32(θ.r_cone),
            Float32(θ.p1_cone),
            Float32(θ.p2_cone),
            Float32(θ.χ_cone),
            Float32(θ.ι_cone),
            Float32(θ.βv_cone),
            Float32(θ.spec_cone),
            Float32(θ.η_cone),
            Float32(θ.r_disk),
            Float32(θ.p1_disk),
            Float32(θ.p2_disk),
            Float32(θ.χ_disk),
            Float32(θ.ι_disk),
            Float32(θ.βv_disk),
            Float32(θ.spec_disk),
            Float32(θ.η_disk),
            Float32(θ.rJ),
            Float32(θ.xtrans),
            Float32(θ.ytrans)
        ), blur
    )
end

function fit(info::ModelInfo, maxevals::Int)
    (;model, name, lower, upper, file, blur) = info
    iobase = joinpath(pwd(), "runs","image_domain", split(file, "/")[end], name)
    mkpath(iobase)
    println("Saving output at :" * iobase)

    inimg = Float32.(VIDA.load_image(file))
    dx, dy = Float32.(rad2μas.([pixelsizes(inimg)...]))
    g = axisdims(inimg)
    vals = zeros(Float32, size(inimg))

    if iszero(blur)
        vals = Float32.(parent(inimg))
    else
        σ_px = blur/(2*sqrt(2*log(2)))/dx
        σ_py = blur/(2*sqrt(2*log(2)))/dy
        
        # Now I need to pick my kernel size. I am going out to 5σ for the
        # gaussian kernel. I have to add one for the convolution to play nice
        nkern = Int(floor(σ_px)*10 + 1)
        # Note I tried to use IIRGaussian but it wasn't accurate enough for us.
        vals = Float32.(imfilter(parent(inimg),
            Kernel.gaussian((σ_py, σ_px),(nkern,nkern)),
            Fill(0.0, inimg),
            Algorithm.FFT()
        ))
    end
    #vals = imfilter(parent(inimg),Kernel.gaussian(blur/(2*√(dx*dy*2log(2f0)))))
    g32 = RectiGrid((X=map(Float32, rad2μas(g.X)), Y=map(Float32, rad2μas(g.Y))))
    img = IntensityMap(vals, g32) #|>MtlArray

    bh = VIDA.NxCorr(img)
    prob = VIDAProblem(bh, model, lower, upper);
    f, t, (lb, ub) = VIDA.build_opt(prob, true)
    fnew(θ) = Float64(f(Float32.(θ)))
    ndim = length(lb)
    resbb =  bboptimize(fnew; NumDimensions=ndim,
                        Method=:adaptive_de_rand_1_bin,
                        MaxFuncEvals=maxevals, TraceMode=:verbose,
                        #PopulationSize=optimizer.popsize,
                        SearchRange=collect(zip(lb, ub)))
                
    xopt = VIDA.transform(t, best_candidate(resbb))
    optfilt =  dual_cone_blur_temp(xopt)
    
    fig = triptic(img, optfilt)

    newkeys = fieldnames(typeof(xopt))
    opt = getfield.(Ref(xopt), newkeys)
    intmap = intensitymap(JuKeBOX(xopt...), μas2rad(160), μas2rad(160), 400, 400) 
    Comrade.save(iobase*"/best.fits", intmap)

    fileout = open(iobase*"/best_nxcorr.txt", "w")

    divmin = divergence(bh, optfilt)
    write(fileout, "nxcorr = "*string(exp(-divmin))*"\n")
    write(fileout, "upper = "*string(upper)*"\n")
    write(fileout, "lower = "*string(lower)*"\n")
    write(fileout, "best_fit = "*string(NamedTuple{newkeys}(opt)))
    close(fileout)

    Plots.savefig(fig, joinpath(iobase,"tripic.png"))

end