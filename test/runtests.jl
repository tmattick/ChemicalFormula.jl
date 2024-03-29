using ChemicalFormula
using Test
using Unitful

@testset "Constructor" begin
    water = Formula("H2O", "water")
    @test water.charge == 0
    @test water.composition == Dict(:H => 2, :O => 1)
    @test water.formula == "H2O"
    @test water.name == "water"

    sulfate = Formula("SO4", -2)
    @test sulfate.charge == -2
    @test sulfate.composition == Dict(:S => 1, :O => 4)
    @test sulfate.formula == "SO4"
    @test isnothing(sulfate.name)

    cyanide = Formula("Fe(CN)6*5H2O")
    @test cyanide.charge == 0
    @test cyanide.composition == Dict(:Fe => 1, :C => 6, :N => 6, :H => 10, :O => 5)
    @test cyanide.formula == "Fe(CN)6*5H2O"
    @test isnothing(cyanide.name)

    nitrate = Formula("NO₃⁻")
    @test nitrate.charge == -1
    @test nitrate.composition == Dict(:N => 1, :O => 3)
    @test nitrate.formula == "NO₃⁻"
    @test isnothing(nitrate.name)

    phosphate = Formula("PO4³⁻", "phosphate")
    @test phosphate.charge == -3
    @test phosphate.composition == Dict(:P => 1, :O => 4)
    @test phosphate.formula == "PO4³⁻"
    @test phosphate.name == "phosphate"

    korund = Formula("Al₂O₃")
    @test korund.charge == 0
    @test korund.composition == Dict(:Al => 2, :O => 3)
    @test korund.formula == "Al₂O₃"
    @test isnothing(korund.name)

    calcium = Formula("Ca²⁺")
    @test calcium.charge == 2
    @test calcium.composition == Dict(:Ca => 1)
    @test calcium.formula == "Ca²⁺"
    @test isnothing(calcium.name)

    potassium = Formula("K⁺")
    @test potassium.charge == 1
    @test potassium.composition == Dict(:K => 1)
    @test potassium.formula == "K⁺"
    @test isnothing(potassium.name)
end

@testset "equality" begin
    @test Formula("H2O") == Formula("H2O")
    @test Formula("NaCl") == Formula("ClNa")
    @test Formula("CCl4", "tetrachloromethane") == Formula("CCl4", "carbon tetrachloride")
    @test Formula("SO4", -2) == Formula("SO4", -2)
    @test Formula("Fe") != Formula("Ni")
    @test Formula("H2O") != Formula("OH")
    @test Formula("Cs") != Formula("Cs", 1)
end

@testset "sumformula" begin
    @test sumformula(Formula("C2H4")) in Set(["C2H4", "H4C2"])
    @test sumformula(Formula("Fe(CN)6", -3)) in
          Set(["FeC6N6", "FeN6C6", "C6FeN6", "C6N6Fe", "N6FeC6", "N6C6Fe"])
    @test sumformula(Formula("Na*6H2O")) in
          Set(["NaH12O6", "NaO6H12", "H12NaO6", "H12O6Na", "O6NaH12", "O6H12Na"])
end

@testset "hillformula" begin
    @test hillformula(Formula("C2H4")) == "C2H4"
    @test hillformula(Formula("H4C2")) == "C2H4"
    @test hillformula(Formula("C2H5OH")) == "C2H6O"
    @test hillformula(Formula("H3O6C6F3")) == "C6H3F3O6"
    @test hillformula(Formula("SO4")) == "O4S"
end

@testset "unicode" begin
    @test unicode(Formula("C2H4")) == "C₂H₄"
    @test unicode(Formula("C2H4"), "sum") in Set(["C₂H₄", "H₄C₂"])
    @test unicode(Formula("Fe(CN)6", -3), false) == "Fe(CN)₆"
    @test unicode(Formula("Fe(CN)6", -3)) == "Fe(CN)₆³⁻"
    @test unicode(Formula("Fe(CN)6", -3), false, "sum") in
          Set(["FeC₆N₆", "FeN₆C₆", "C₆FeN₆", "C₆N₆Fe", "N₆FeC₆", "N₆C₆Fe"])
    @test unicode(Formula("Fe(CN)6", -3), "sum") in
          Set(["FeC₆N₆³⁻", "FeN₆C₆³⁻", "C₆FeN₆³⁻", "C₆N₆Fe³⁻", "N₆FeC₆³⁻", "N₆C₆Fe³⁻"])
    @test unicode(Formula("H4C2"), "hill") == "C₂H₄"
    @test unicode(Formula("H3O6C6F3", -3), "hill") == "C₆H₃F₃O₆³⁻"
end

@testset "latex" begin
    @test latex(Formula("C2H4")) == raw"\ce{C2H4}"
    @test latex(Formula("C2H4"), "sum") in Set([raw"\ce{C2H4}", raw"\ce{H4C2}"])
    @test latex(Formula("Fe(CN)6", -3), false) == raw"\ce{Fe(CN)6}"
    @test latex(Formula("Fe(CN)6", -3)) == raw"\ce{Fe(CN)6^{3-}}"
    @test latex(Formula("Fe(CN)6", -3), false, "sum") in Set([
        raw"\ce{FeC6N6}",
        raw"\ce{FeN6C6}",
        raw"\ce{C6FeN6}",
        raw"\ce{C6N6Fe}",
        raw"\ce{N6FeC6}",
        raw"\ce{N6C6Fe}",
    ])
    @test latex(Formula("Fe(CN)6", -3), "sum") in Set([
        raw"\ce{FeC6N6^{3-}}",
        raw"\ce{FeN6C6^{3-}}",
        raw"\ce{C6FeN6^{3-}}",
        raw"\ce{C6N6Fe^{3-}}",
        raw"\ce{N6FeC6^{3-}}",
        raw"\ce{N6C6Fe^{3-}}",
    ])
    @test latex(Formula("H4C2"), "hill") == raw"\ce{C2H4}"
    @test latex(Formula("H3O6C6F3", -3), "hill") == raw"\ce{C6H3F3O6^{3-}}"
end

@testset "formulamass" begin
    @test formulamass(Formula("Fe")) ≈ 55.845u"u"
    @test formulamass(Formula("PO4", -3)) ≈ 94.9697619985u"u"
    @test formulamass(Formula("C2H6", "ethane")) ≈ 30.07u"u"
end

@testset "formulaweight" begin
    @test formulaweight(Formula("U")) ≈ 238.028913u"g/mol"
    @test formulaweight(Formula("NH4", 1)) ≈ 18.039u"g/mol"
    @test formulaweight(Formula("CO2H2", "formic acid")) ≈ 46.025u"g/mol"
end

@testset "massfractions" begin
    @test massfractions(Formula("C")) == Dict(:C => 1.0)
    @test massfractions(Formula("NaCl")) ==
          Dict(:Na => 0.3933925401014177, :Cl => 0.6066074598985822)
    @test massfractions(Formula("C2H4")) ==
          Dict(:C => 0.8562771797248164, :H => 0.14372282027518357)
end

@testset "isradioactive" begin
    @test isradioactive(Formula("Ra")) == true
    @test isradioactive(Formula("Mn")) == false
    @test isradioactive(Formula("TbUC2")) == true
    @test isradioactive(Formula("Fe(CN)6")) == false
end

@testset "ischarged" begin
    @test ischarged(Formula("H2O")) == false
    @test ischarged(Formula("Na", 1)) == true
    @test ischarged(Formula("S", -2)) == true
end

@testset "textcharge" begin
    @test textcharge(Formula("CH4")) == ""
    @test textcharge(Formula("Li", 1)) == "+"
    @test textcharge(Formula("F", -1)) == "-"
    @test textcharge(Formula("Ca", 2)) == "2+"
    @test textcharge(Formula("PO4", -3)) == "3-"
end
