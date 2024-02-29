struct MarginMakieHist <: PairPlots.VizTypeDiag
    kwargs
    MarginMakieHist(;kwargs...) = new(kwargs)
end

function PairPlots.diagplot(ax::Makie.Axis, viz::MarginMakieHist, series::PairPlots.AbstractSeries, colname)

    cn = PairPlots.columnnames(series)
    if colname ∉ cn
        return
    end
    dat = getproperty(series.table, colname)

    bins = get(series.kwargs, :bins, 32)
    bins = get(viz.kwargs, :bins, bins)
    
    Makie.hist!(ax, dat; 
    series.kwargs..., 
    viz.kwargs..., 
    bins=bins)
    Makie.ylims!(ax,low=0)
end

