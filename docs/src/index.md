# ChemicalFormula.jl Documentation

ChemicalFormula.jl is a package for the simple representation of chemical compounds as 
formulas. The package aims to be a lightweight solution for handling various chemical 
formulas while providing general information and formatting capabilities. It is strongly 
inspired by the Python package [ChemFormula](https://github.com/molshape/ChemFormula).

The focus of ChemicalFormula.jl is on the parsing capabilities. These include the handling 
of Parantheses and coordinating molecules marked with "*". From these you can generate 
generate formatted strings (text, Unicode and LaTeX) in various styles.

In addition, ChemicalFormula.jl offers the calculation of formula weights and thus enables 
further stochiometric calculations. The atomic weights correspond to the specifications of 
the IUPAC Commission for Isotope abundances and atomic weights. The values are taken from 
https://iupac.qmul.ac.uk/AtWt/.

## Installation

ChemicalFormula can be installed using the package manager. Enter the Pkg REPL mode by 
typing `]` in the REPL and run
```
pkg> add ChemicalFormula
```

## Examples

```julia-repl
julia> using ChemicalFormula
julia> water = Formula("H2O", "water")
Formula("H2O", Dict{Symbol, Int32}(:H => 2, :O => 1), 0, "water")
julia> formulaweight(water)
18.015 g mol⁻¹
julia> perfluorotrimesate = Formula("H3O6C6F3", -3)
Formula("H3O6C6F3", Dict{Symbol, Int32}(:F => 3, :H => 3, :O => 6, :C => 6), -3, nothing)
julia> unicode(perfluorotrimesate, "hill")
"C₆H₃F₃O₆³⁻"
julia> charged(perfluorotrimesate)
true
julia> cyanide = Formula("K4Fe(CN)6")
Formula("K4Fe(CN)6", Dict{Symbol, Int32}(:N => 6, :Fe => 1, :K => 4, :C => 6), 0, nothing)
julia> latex(cyanide, "sum")
"\\ce{N6FeK4C6}"
julia> charged(cyanide)
false
julia> uranyl = Formula("UO2F2*H2O")
Formula("UO2F2*H2O", Dict{Symbol, Int32}(:U => 1, :F => 2, :H => 2, :O => 3), 0, nothing)
julia> radioactive(uranyl)
true
```
