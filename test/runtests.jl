using ChemicalFormula
using Test
using PeriodicTable, Unitful

@testset "Constructor" begin
    water = Formula("H2O", "water")
    @test water.charge == 0
    @test water.composition == Dict(elements[:H] => 2, elements[:O] => 1)
    @test water.formula == "H2O"
    @test water.name == "water"

    sulfate = Formula("SO4", -2)
    @test sulfate.charge == -2
    @test sulfate.composition == Dict(elements[:S] => 1, elements[:O] => 4)
    @test sulfate.formula == "SO4"
    @test isnothing(sulfate.name)

    cyanide = Formula("Fe(CN)6*5H2O")
    @test cyanide.charge == 0
    @test cyanide.composition == Dict(
        elements[:Fe] => 1,
        elements[:C] => 6,
        elements[:N] => 6,
        elements[:H] => 10,
        elements[:O] => 5,
    )
    @test cyanide.formula == "Fe(CN)6*5H2O"
    @test isnothing(cyanide.name)
end

@testset "formulamass" begin
    @test formulamass(Formula("Fe")) ≈ 55.8452u"u"
    @test formulamass(Formula("PO4", -3)) ≈ 94.9697619985u"u"
    @test formulamass(Formula("C2H6", "ethane")) ≈ 30.07u"u"
end

@testset "formulaweight" begin
    @test formulaweight(Formula("U")) ≈ 238.028913u"g/mol"
    @test formulaweight(Formula("NH4", 1)) ≈ 18.039u"g/mol"
    @test formulaweight(Formula("CO2H2", "formic acid")) ≈ 46.025u"g/mol"
end

@testset "massfractions" begin
    @test massfractions(Formula("C")) == Dict(elements[:C] => 1.0)
    @test massfractions(Formula("NaCl")) ==
          Dict(elements[:Na] => 0.39339254012217784, elements[:Cl] => 0.6066074598778223)
    @test massfractions(Formula("C2H4")) ==
          Dict(elements[:C] => 0.8562771797248164, elements[:H] => 0.14372282027518357)
end

@testset "radioactive" begin
    @test radioactive(Formula("Ra")) == true
    @test radioactive(Formula("Mn")) == false
    @test radioactive(Formula("TbUC2")) == true
    @test radioactive(Formula("Fe(CN)6")) == false
end

@testset "charged" begin
    @test charged(Formula("H2O")) == false
    @test charged(Formula("Na", 1)) == true
    @test charged(Formula("S", -2)) == true
end

@testset "textcharge" begin
    @test textcharge(Formula("CH4")) == ""
    @test textcharge(Formula("Li", 1)) == "+"
    @test textcharge(Formula("F", -1)) == "-"
    @test textcharge(Formula("Ca", 2)) == "2+"
    @test textcharge(Formula("PO4", -3)) == "3-"
end
