using GrowthTop
using Documenter
using PUPM

DocMeta.setdocmeta!(GrowthTop, :DocTestSetup, :(using GrowthTop, PUPM); recursive=true)

makedocs(;
    modules=[GrowthTop],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    sitename="GrowthTop.jl",
    format=Documenter.HTML(;
        canonical="https://Aminofa70.github.io/GrowthTop.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Aminofa70/GrowthTop.jl",
    devbranch="main",
)
