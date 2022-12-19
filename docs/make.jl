using Documenter
using ChemicalFormula

makedocs(
    sitename = "ChemicalFormula.jl",
    modules = [ChemicalFormula],
    pages = Any["Home"=>"index.md", "Call Reference"=>"reference.md"],
)

deploydocs(
    repo = "github.com/tmattick/ChemicalFormula.jl.git",
)

