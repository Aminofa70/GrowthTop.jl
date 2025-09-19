using Pkg
Pkg.add(url = "https://github.com/Aminofa70/PUPM.jl", rev= "develop"); 
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()
# cd(@__DIR__)
# Pkg.activate(".")
# push!(LOAD_PATH,"../src/")
using GrowthTop
#@show mm # just to test that the package is loaded
using Documenter

DocMeta.setdocmeta!(GrowthTop, :DocTestSetup, :(using GrowthTop); recursive=true)


makedocs(;
    modules=[GrowthTop.jl],
    authors="Aminofa70 <amin.alibakhshi@upm.es> and contributors",
    repo="https://github.com/Aminofa70/GrowthTop.jl/blob/{commit}{path}#L{line}",
    sitename="GrowthTop.jl",
    # format=Documenter.HTML(;
        # prettyurls=get(ENV, "CI", "false") == "true",
        # canonical="https://olejorik.github.io/PhaseRetrieval.jl",
        # assets=String[],
    # ),
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical = "https://Aminofa70.github.io/GrowthTop.jl",
        assets = ["assets/favicon.ico"],
        highlights = ["yaml"],
    ),
    clean = false,
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="https://Aminofa70.github.io/GrowthTop.jl",
    target = "build",
)

Pkg.rm("PUPM.jl")