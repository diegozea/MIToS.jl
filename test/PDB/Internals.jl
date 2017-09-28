@testset "AtomsData internals" begin

    @test length(PDB._3_letter_aa) == 20

    s = Set{Tuple{String,String}}()
    PDB._add_CTER_O!(s)
    @test length(s) == 20 * 3
    for aa in PDB._3_letter_aa
        @test (aa,"OXT") in s
        @test (aa,"OT2") in s
        @test (aa,"OT1") in s
    end

    d = Dict{Tuple{String,String},Float64}()
    PDB._add_CTER_O!(d, 1.0)
    @test length(d) == 20 * 3
    for aa in PDB._3_letter_aa
        @test d[(aa,"OXT")] == 1.0
        @test d[(aa,"OT2")] == 1.0
        @test d[(aa,"OT1")] == 1.0
    end

    value = Set(String["OXT", "OT2", "OT1"])

    gd = Dict{String,Set{String}}()
    PDB._generate_dict!(gd, d)
    @test length(gd) == 20
    for aa in PDB._3_letter_aa
        @test gd[aa] == value
    end

    gs = Dict{String,Set{String}}()
    PDB._generate_dict!(gs, d)
    @test length(gs) == 20
    for aa in PDB._3_letter_aa
        @test gs[aa] == value
    end

    as = PDB._generate_interaction_keys(d, s, s, s, s)
    @test length(as) == 20
    for aa in PDB._3_letter_aa
        @test as[aa] == value
    end

    @test isa(PDB._interaction_keys, Dict{String,Set{String}})
end
