using BayesianBlackHolePhotogrammetry
using Documenter

DocMeta.setdocmeta!(BayesianBlackHolePhotogrammetry, :DocTestSetup, :(using BayesianBlackHolePhotogrammetry); recursive=true)

makedocs(;
    modules=[BayesianBlackHolePhotogrammetry],
    authors="Dominic <dchang3419@hotmail.com> and contributors",
    sitename="BayesianBlackHolePhotogrammetry.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
