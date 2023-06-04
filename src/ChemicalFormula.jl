module ChemicalFormula

export Formula,
    sumformula,
    hillformula,
    unicode,
    latex,
    formulamass,
    formulaweight,
    massfractions,
    isradioactive,
    ischarged,
    textcharge

import Base.==
using Unitful
include("elements.jl")
using .Elements
include("subsuperscript.jl")
using .SubSuperScript

"""
    Formula(formula[, composition], charge=0, name=nothing)

Represent a chemical `formula` with an electrical `charge` and an optional `name`.
The `composition` is automatically determined from the specified `formula`. Compounds can
be grouped with parentheses, coordinating molecules can be annotated with a *. Elements are
represented by `Symbol`s.
If the `charge` is not specified, it is determined from the charge in superscript 
within the `formula`.

# Examples
```julia-repl
julia> Formula("H2O", "water")
Formula("H2O", Dict{Symbol, Int32}(:H => 2, :O => 1), 0, "water")
julia> Formula("SO4", -2)
Formula("SO4", Dict{Symbol, Int32}(:S => 1, :O => 4), -2, nothing)
julia> Formula("Fe(CN)6*5H2O")
Formula("Fe(CN)6*5H2O", Dict{Symbol, Int32}(:N => 6, :Fe => 1, :H => 10, :O => 5, :C => 6), 0, nothing)
julia> Formula("SO₄²⁻")
Formula("SO₄²⁻", Dict{Symbol, Int32}(:S => 1, :O => 4), -2, nothing)
```
"""
struct Formula
    formula::AbstractString
    composition::Dict{Symbol,Int32}
    charge::Int8
    name::Union{AbstractString,Nothing}
end

Formula(formula::AbstractString, charge::Integer, name::AbstractString) =
    Formula(formula, parseformula(formula), charge, name)

Formula(formula::AbstractString, charge::Integer) =
    Formula(formula, parseformula(formula), charge, nothing)

function Formula(formula::AbstractString, name::AbstractString)
    (composition, charge) = parseformulawithcharge(formula)
    return Formula(formula, composition, charge, name)
end

function Formula(formula::AbstractString)
    (composition, charge) = parseformulawithcharge(formula)
    return Formula(formula, composition, charge, nothing)
end

==(f1::Formula, f2::Formula) = f1.composition == f2.composition && f1.charge == f2.charge

"Give a collapsed sum formula with no special ordering and without charge as unformatted text."
function sumformula(formula::Formula)
    out = ""
    for (element, count) in formula.composition
        out *= string(element) * (count == 1 ? "" : string(count))
    end
    return out
end

"Give a sum formula in Hill notation without charge as unformatted text."
function hillformula(formula::Formula)
    out = ""
    elementssorted = sort(collect(keys(formula.composition)))
    if :C in elementssorted
        out *= "C" * (formula.composition[:C] == 1 ? "" : string(formula.composition[:C]))
        deleteat!(elementssorted, findall(x -> x === :C, elementssorted))
        if :H in elementssorted
            out *=
                "H" * (formula.composition[:H] == 1 ? "" : string(formula.composition[:H]))
            deleteat!(elementssorted, findall(x -> x === :H, elementssorted))
        end
    end
    for element in elementssorted
        out *=
            string(element) *
            (formula.composition[element] == 1 ? "" : string(formula.composition[element]))
    end
    return out
end

const displayforms = Dict(
    "formula" => f -> f.formula,
    "hill" => hillformula,
    "hillformula" => hillformula,
    "sum" => sumformula,
    "sumformula" => sumformula,
)

"""
    unicode(formula, includecharge=true, form="formula")

Give a formula written in unicode characters. `includecharge` determines whether the 
electronic charge should be given, defaults to `true`. `form` has to be one of `"formula"`
(the `formula` that was given to construct the `Formula` object), `"hill"` or 
`"hillformula"` for Hill notation or `"sum"` or `"sumformula"` for a collapsed sum 
formula, defaults to `"formula"`.
"""
function unicode(formula::Formula, includecharge::Bool, form::AbstractString)
    sub = Dict(
        '0' => '₀',
        '1' => '₁',
        '2' => '₂',
        '3' => '₃',
        '4' => '₄',
        '5' => '₅',
        '6' => '₆',
        '7' => '₇',
        '8' => '₈',
        '9' => '₉',
    )
    super = Dict(
        '0' => '⁰',
        '1' => '¹',
        '2' => '²',
        '3' => '³',
        '4' => '⁴',
        '5' => '⁵',
        '6' => '⁶',
        '7' => '⁷',
        '8' => '⁸',
        '9' => '⁹',
        '+' => '⁺',
        '-' => '⁻',
    )

    out = displayforms[form](formula)
    out = replace(out, sub...)

    if includecharge
        charge = textcharge(formula)
        charge = replace(charge, super...)
        out *= charge
    end

    return out
end

unicode(formula::Formula, includecharge::Bool) = unicode(formula, includecharge, "formula")

unicode(formula::Formula, form::AbstractString) = unicode(formula, true, form)

unicode(formula::Formula) = unicode(formula, true, "formula")

"""
    latex(formula, includecharge=true, form="formula")

Give a formula written in LaTeX, using the mhchem package. `includecharge` determines 
whether the electronic charge should be given, defaults to `true`. `form` has to be one of 
`"formula"` (the `formula` that was given to construct the `Formula` object), `"hill"` or 
`"hillformula"` for Hill notation or `"sum"` or `"sumformula"` for a collapsed sum 
formula, defaults to `"formula"`.
"""
function latex(formula::Formula, includecharge::Bool, form::AbstractString)
    out = "\\ce{" * displayforms[form](formula)

    if includecharge
        charge = textcharge(formula)
        out *= charge == "" ? "" : "^{$charge}"
    end

    return out * "}"
end

latex(formula::Formula, includecharge::Bool) = latex(formula, includecharge, "formula")

latex(formula::Formula, form::AbstractString) = latex(formula, true, form)

latex(formula::Formula) = latex(formula, true, "formula")

"Return the mass of one `formula` without specified unit."
function mass_wo_unit(formula::Formula)
    mass = 0.0
    for (element, count) in formula.composition
        mass += count * atomicmass[element]
    end
    return mass
end

"Return the mass of one `formula` unit in unified atomic mass units ``u``."
formulamass(formula::Formula) = mass_wo_unit(formula) * 1u"u"

"Return the molar mass of the `formula` in ``g mol⁻¹``."
formulaweight(formula::Formula) = mass_wo_unit(formula) * 1u"g/mol"

"Return the mass fraction for each element of the `formula` in a `element => fraction` `Dict`."
function massfractions(formula::Formula)
    totalmass = mass_wo_unit(formula)
    fractions = Dict{Symbol,Float64}()
    for (element, count) in formula.composition
        fraction = (count * atomicmass[element]) / totalmass
        fractions[element] = fraction
    end
    return fractions
end

"Return whether the `formula` is radioactive, i.e. contains radioactive elements."
function isradioactive(formula::Formula)
    for element in keys(formula.composition)
        elementnumber = atomicnumbers[element]
        if elementnumber >= 84 || elementnumber == 61 || elementnumber == 43
            return true
        end
    end
    return false
end

"Return whether the `formula` carries an electrical charge."
ischarged(formula::Formula) = formula.charge != 0

"Return formatted `String` of the `formula` `charge` (like 4+, 3-, +, ...)."
function textcharge(formula::Formula)
    if formula.charge == 0
        return ""
    elseif formula.charge == 1
        return "+"
    elseif formula.charge == -1
        return "-"
    elseif formula.charge > 0
        return "$(formula.charge)+"
    else
        return "$(-formula.charge)-"
    end
end

"Parse a formula `String` into a `Dict` with corresponding element, count pairs."
function parseformula(s::AbstractString)
    s = reduceformula(s)
    composition = Dict{Symbol,Int32}()
    for m in eachmatch(r"(?<element>[A-Z][a-z]?)(?<count>\d*)", s)
        element = Symbol(m["element"])
        num = m["count"] == "" ? one(Int32) : parse(Int32, m["count"])
        if element in keys(composition)
            composition[element] += num
        else
            composition[element] = num
        end
    end
    return composition
end

"""
Parse a formula `String` into a `Dict` with corresponding element, count pairs while parsing
the `charge` from superscript. Defaults to `parsefomula(s)` if no superscript is present.
"""
function parseformulawithcharge(s::AbstractString)
    s = reduceformula(s)

    charge = zero(Int8)
    i = 1
    if issuperscript(s[end])
        while i <= length(s) && !issuperscript(s[i])
            i += 1
        end
        charge = parsecharge(s[i:end])
    else
        return parseformula(s), zero(Int8)
    end

    s = s[1:i-1]
    composition = Dict{Symbol,Int32}()
    for m in eachmatch(r"(?<element>[A-Z][a-z]?)(?<count>\d*)", s)
        element = Symbol(m["element"])
        num = m["count"] == "" ? one(Int32) : parse(Int32, m["count"])
        if element in keys(composition)
            composition[element] += num
        else
            composition[element] = num
        end
    end
    return composition, charge
end

"Reduces a formula `String` by removing whitespace, subscripts, asterisks and parentheses."
function reduceformula(s::AbstractString)
    s = replace(s, r"\s" => "")
    s = subscripttonormal(s)
    s = removestar(s)
    s = removebrackets(s)
    return s
end

"Remove '*' `char`s in a formula `String` by rewriting the `String` as a single molecule formula."
function removestar(s::AbstractString)
    if !occursin('*', s)
        return s
    end

    (beforestar, afterstar) = rsplit(s, "*", limit = 2)
    afterunit = match(r"(?<multiplier>\d*)(?<formula>[A-Za-z0-9]+)", afterstar)
    aftercomponents =
        eachmatch(r"(?<element>[A-Z][a-z]?)(?<count>\d*)", afterunit["formula"])
    aftermultiplier =
        afterunit["multiplier"] == "" ? 1 : parse(Int, afterunit["multiplier"])
    afterformula = ""
    for comp in aftercomponents
        frequency = comp["count"] == "" ? 1 : parse(Int, comp["count"])
        afterformula *= comp["element"] * string(frequency * aftermultiplier)
    end
    return removestar(beforestar * afterformula)
end

"Remove parentheses in a formula by rewriting it without grouped parts of a compound."
function removebrackets(s::AbstractString)
    if !occursin('(', s)
        return s
    end

    innerunit = match(r"\((?<inbrackets>[A-Za-z0-9]+)\)(?<bracketmultiplier>\d*)", s) # inside brackets with no brackets inside
    beforebrackets = s[1:innerunit.offsets[1]-2]
    afterbrackets = s[innerunit.offsets[2]+length(innerunit.captures[2]):end]
    innercomponents =
        eachmatch(r"(?<element>[A-Z][a-z]?)(?<count>\d*)", innerunit["inbrackets"])
    innermultiplier =
        innerunit["bracketmultiplier"] == "" ? 1 :
        parse(Int, innerunit["bracketmultiplier"])
    innerformula = ""
    for comp in innercomponents
        frequency = comp["count"] == "" ? 1 : parse(Int, comp["count"])
        innerformula *= comp["element"] * string(frequency * innermultiplier)
    end
    return removebrackets(beforebrackets * innerformula * afterbrackets)
end

"Parse the charge of a formula given in superscript."
function parsecharge(s::AbstractString)
    s = superscripttonormal(s)
    if s == "+"
        return one(Int8)
    elseif s == "-"
        return -one(Int8)
    else
        numcharge = parse(Int8, s[1:end-1])
        if s[end] == '+'
            return numcharge
        else
            return -numcharge
        end
    end
end

end
