using SPM
using Documenter

DocMeta.setdocmeta!(SPM, :DocTestSetup, :(using SPM); recursive=true)

makedocs(;
  modules=[SPM],
  authors="SPM team and contributors",
  repo="https://github.com/spm/SPM.jl/blob/{commit}{path}#{line}",
  sitename="SPM.jl",
  doctest=true,
  format=Documenter.HTML(;
    prettyurls=get(ENV, "Tests", "false") == "true",
    canonical="https://spm.github.io/SPM.jl",
    edit_link="main",
    assets=String[]
  ),
  pages=[
    "Introduction" => "index.md",
    "Installation" => "installation.md",
    "Example" => Any[
            "Coregistration" => "coregistration.md",],
    "API" =>"API.md",
  ]
)

deploydocs(;
  repo="github.com/spm/SPM.jl",
  devbranch="main"
)
