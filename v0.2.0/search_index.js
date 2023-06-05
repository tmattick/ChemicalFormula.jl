var documenterSearchIndex = {"docs":
[{"location":"reference/#Call-Reference","page":"Call Reference","title":"Call Reference","text":"","category":"section"},{"location":"reference/#Constructor","page":"Call Reference","title":"Constructor","text":"","category":"section"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"Formula","category":"page"},{"location":"reference/#ChemicalFormula.Formula","page":"Call Reference","title":"ChemicalFormula.Formula","text":"Formula(formula[, composition], charge=0, name=nothing)\n\nRepresent a chemical formula with an electrical charge and an optional name. The composition is automatically determined from the specified formula. Compounds can be grouped with parentheses, coordinating molecules can be annotated with a *. Elements are represented by Symbols.\n\nExamples\n\njulia> Formula(\"H2O\", \"water\")\nFormula(\"H2O\", Dict{Symbol, Int32}(:H => 2, :O => 1), 0, \"water\")\njulia> Formula(\"SO4\", -2)\nFormula(\"SO4\", Dict{Symbol, Int32}(:S => 1, :O => 4), -2, nothing)\njulia> Formula(\"Fe(CN)6*5H2O\")\nFormula(\"Fe(CN)6*5H2O\", Dict{Symbol, Int32}(:N => 6, :Fe => 1, :H => 10, :O => 5, :C => 6), 0, nothing)\n\n\n\n\n\n","category":"type"},{"location":"reference/#String-formatting","page":"Call Reference","title":"String formatting","text":"","category":"section"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"unicode","category":"page"},{"location":"reference/#ChemicalFormula.unicode","page":"Call Reference","title":"ChemicalFormula.unicode","text":"unicode(formula, includecharge=true, form=\"formula\")\n\nGive a formula written in unicode characters. includecharge determines whether the  electronic charge should be given, defaults to true. form has to be one of \"formula\" (the formula that was given to construct the Formula object), \"hill\" or  \"hillformula\" for Hill notation or \"sum\" or \"sumformula\" for a collapsed sum  formula, defaults to \"formula\".\n\n\n\n\n\n","category":"function"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"latex","category":"page"},{"location":"reference/#ChemicalFormula.latex","page":"Call Reference","title":"ChemicalFormula.latex","text":"latex(formula, includecharge=true, form=\"formula\")\n\nGive a formula written in LaTeX, using the mhchem package. includecharge determines  whether the electronic charge should be given, defaults to true. form has to be one of  \"formula\" (the formula that was given to construct the Formula object), \"hill\" or  \"hillformula\" for Hill notation or \"sum\" or \"sumformula\" for a collapsed sum  formula, defaults to \"formula\".\n\n\n\n\n\n","category":"function"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"sumformula","category":"page"},{"location":"reference/#ChemicalFormula.sumformula","page":"Call Reference","title":"ChemicalFormula.sumformula","text":"Give a collapsed sum formula with no special ordering and without charge as unformatted text.\n\n\n\n\n\n","category":"function"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"hillformula","category":"page"},{"location":"reference/#ChemicalFormula.hillformula","page":"Call Reference","title":"ChemicalFormula.hillformula","text":"Give a sum formula in Hill notation without charge as unformatted text.\n\n\n\n\n\n","category":"function"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"textcharge","category":"page"},{"location":"reference/#ChemicalFormula.textcharge","page":"Call Reference","title":"ChemicalFormula.textcharge","text":"Return formatted String of the formula charge (like 4+, 3-, +, ...).\n\n\n\n\n\n","category":"function"},{"location":"reference/#Mass-and-weight-related-functions","page":"Call Reference","title":"Mass and weight related functions","text":"","category":"section"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"formulamass","category":"page"},{"location":"reference/#ChemicalFormula.formulamass","page":"Call Reference","title":"ChemicalFormula.formulamass","text":"Return the mass of one formula unit in unified atomic mass units u.\n\n\n\n\n\n","category":"function"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"formulaweight","category":"page"},{"location":"reference/#ChemicalFormula.formulaweight","page":"Call Reference","title":"ChemicalFormula.formulaweight","text":"Return the molar mass of the formula in g mol¹.\n\n\n\n\n\n","category":"function"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"massfractions","category":"page"},{"location":"reference/#ChemicalFormula.massfractions","page":"Call Reference","title":"ChemicalFormula.massfractions","text":"Return the mass fraction for each element of the formula in a element => fraction Dict.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Utility-functions","page":"Call Reference","title":"Utility functions","text":"","category":"section"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"ischarged","category":"page"},{"location":"reference/#ChemicalFormula.ischarged","page":"Call Reference","title":"ChemicalFormula.ischarged","text":"Return whether the formula carries an electrical charge.\n\n\n\n\n\n","category":"function"},{"location":"reference/","page":"Call Reference","title":"Call Reference","text":"isradioactive","category":"page"},{"location":"reference/#ChemicalFormula.isradioactive","page":"Call Reference","title":"ChemicalFormula.isradioactive","text":"Return whether the formula is radioactive, i.e. contains radioactive elements.\n\n\n\n\n\n","category":"function"},{"location":"#ChemicalFormula.jl-Documentation","page":"Home","title":"ChemicalFormula.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ChemicalFormula.jl is a package for the simple representation of chemical compounds as  formulas. The package aims to be a lightweight solution for handling various chemical  formulas while providing general information and formatting capabilities. It is strongly  inspired by the Python package ChemFormula.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The focus of ChemicalFormula.jl is on the parsing capabilities. These include the handling  of Parantheses and coordinating molecules marked with \"*\". From these you can generate  generate formatted strings (text, Unicode and LaTeX) in various styles.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In addition, ChemicalFormula.jl offers the calculation of formula weights and thus enables  further stochiometric calculations. The atomic weights correspond to the specifications of  the IUPAC Commission for Isotope abundances and atomic weights. The values are taken from  https://iupac.qmul.ac.uk/AtWt/.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ChemicalFormula can be installed using the package manager. Enter the Pkg REPL mode by  typing ] in the REPL and run","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add ChemicalFormula","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using ChemicalFormula\njulia> water = Formula(\"H2O\", \"water\")\nFormula(\"H2O\", Dict{Symbol, Int32}(:H => 2, :O => 1), 0, \"water\")\njulia> formulaweight(water)\n18.015 g mol⁻¹\njulia> perfluorotrimesate = Formula(\"H3O6C6F3\", -3)\nFormula(\"H3O6C6F3\", Dict{Symbol, Int32}(:F => 3, :H => 3, :O => 6, :C => 6), -3, nothing)\njulia> unicode(perfluorotrimesate, \"hill\")\n\"C₆H₃F₃O₆³⁻\"\njulia> ischarged(perfluorotrimesate)\ntrue\njulia> cyanide = Formula(\"K4Fe(CN)6\")\nFormula(\"K4Fe(CN)6\", Dict{Symbol, Int32}(:N => 6, :Fe => 1, :K => 4, :C => 6), 0, nothing)\njulia> latex(cyanide, \"sum\")\n\"\\\\ce{N6FeK4C6}\"\njulia> ischarged(cyanide)\nfalse\njulia> uranyl = Formula(\"UO2F2*H2O\")\nFormula(\"UO2F2*H2O\", Dict{Symbol, Int32}(:U => 1, :F => 2, :H => 2, :O => 3), 0, nothing)\njulia> isradioactive(uranyl)\ntrue","category":"page"}]
}