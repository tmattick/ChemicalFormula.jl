using ChemicalFormula
using Test
using PeriodicTable

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
