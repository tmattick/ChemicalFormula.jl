module ChemicalFormula

export Formula

import Base.==
using PeriodicTable, Unitful

"""
    Formula(formula[, composition], charge=0, name=nothing)

Represent a chemical `formula` with an electrical `charge` and an optional `name`.
The `composition` is automatically determined from the specified `formula`. Compounds can
be grouped with parentheses, coordinating molecules can be annotated with a *. Elements are
represented by the `PeriodicTable.Element` type.

# Examples
```julia-repl
julia> Formula("H2O", "water")
Formula("H2O", Dict{PeriodicTable.Element, UInt32}(Element(Oxygen) => 0x00000001, Element(Hydrogen) => 0x00000002), 0, "water")
julia> Formula("SO4", -2)
Formula("SO4", Dict{PeriodicTable.Element, UInt32}(Element(Sulfur) => 0x00000001, Element(Oxygen) => 0x00000004), -2, nothing)
julia> Formula("Fe(CN)6*5H2O")
Formula("Fe(CN)6*5H2O", Dict{PeriodicTable.Element, UInt32}(Element(Carbon) => 0x00000006, Element(Nitrogen) => 0x00000006, Element(Iron) => 0x00000001, Element(Oxygen) => 0x00000005, Element(Hydrogen) => 0x0000000a), 0, nothing)
```
"""
struct Formula
    formula::AbstractString
    composition::Dict{Element,UInt32}
    charge::Int8
    name::Union{AbstractString,Nothing}
end

Formula(formula::AbstractString, charge::Integer, name::AbstractString) =
    Formula(formula, parseformula(formula), charge, name)

Formula(formula::AbstractString, charge::Integer) =
    Formula(formula, parseformula(formula), charge, nothing)

Formula(formula::AbstractString, name::AbstractString) =
    Formula(formula, parseformula(formula), 0, name)

Formula(formula::AbstractString) = Formula(formula, parseformula(formula), 0, nothing)

==(f1::Formula, f2::Formula) = f1.composition == f2.composition && f1.charge == f2.charge

"Return the mass of one `formula` unit in unified atomic mass units ``u``."
function formulamass(formula::Formula)
    mass = 0.0u"u"
    for (element, count) in formula.composition
        mass += count * element.atomic_mass
    end
    return mass
end

"Return the molar mass of the `formula` in ``g mol⁻¹``."
formulaweight(formula::Formula) = formulamass(formula) * 1u"g/(mol*u)"

"Return the mass fraction for each element of the `formula` in a `element => fraction` Dict."
function massfractions(formula::Formula)
    totalmass = formulamass(formula)
    fractions = Dict{Element,Float64}()
    for (element, count) in formula.composition
        fraction = (count * element.atomic_mass) / totalmass
        fractions[element] = fraction
    end
    return fractions
end

"Return whether the `formula` is radioactive, i.e. contains radioactive elements."
function radioactive(formula::Formula)
    for element in keys(formula.composition)
        if element.number >= 84 || element.number == 61 || element.number == 43
            return true
        end
    end
    return false
end

"Return whether the `formula` carries an electrical charge."
charged(formula::Formula) = formula.charge != 0

"Return formatted String of the `formula` `charge` (like 4+, 3-, +, ...)."
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

"Parse a formula String into a dictionary with corresponding element, count pairs."
function parseformula(s::AbstractString)
    s = replace(s, " " => "")
    s = removestar(s)
    s = removebrackets(s)
    composition = Dict{Element,UInt32}()
    for m in eachmatch(r"(?<element>[A-Z][a-z]?)(?<count>\d*)", s)
        element = elements[Symbol(m["element"])]
        num = m["count"] == "" ? one(UInt32) : parse(UInt32, m["count"])
        if element in keys(composition)
            composition[element] += num
        else
            composition[element] = num
        end
    end
    return composition
end

"Remove * characters in a formula string by rewriting the string as a single molecule formula."
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

end
