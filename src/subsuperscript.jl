module SubSuperScript

export issuperscript, subscripttonormal, superscripttonormal

const superscripts = Dict{Char, Char}(
    '⁰' => '0',
    '¹' => '1',
    '²' => '2',
    '³' => '3',
    '⁴' => '4',
    '⁵' => '5',
    '⁶' => '6',
    '⁷' => '7',
    '⁸' => '8',
    '⁹' => '9',
    '⁺' => '+',
    '⁻' => '-',
)

const subscripts = Dict{Char, Char}(
    '₀' => '0',
    '₁' => '1',
    '₂' => '2',
    '₃' => '3',
    '₄' => '4',
    '₅' => '5',
    '₆' => '6',
    '₇' => '7',
    '₈' => '8',
    '₉' => '9',
)

"Return whether `c` is a numeric superscript or ⁺/⁻."
issuperscript(c::Char) = c in keys(superscripts)

"Return whether `c` is a numeric subscript."
issubscript(c::Char) = c in keys(subscripts)

"Convert all numeric superscripts or ⁺/⁻ in `s` to normal line."
superscripttonormal(s::AbstractString) = replace(s, superscripts...)

"Convert all numeric subscripts in `s` to normal line."
subscripttonormal(s::AbstractString) = replace(s, subscripts...)

end