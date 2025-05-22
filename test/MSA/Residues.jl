using MIToS.MSA # Imports Residue, GAP, XAA, etc.
using Test

@testset "Residue Show/Print" begin
    @testset "Base.show(io::IO, res::Residue)" begin
        buf = IOBuffer()
        show(buf, Residue('A'))
        @test String(take!(buf)) == "A"
        
        # Reset buffer for next test
        # buf = IOBuffer() # Not needed as take! clears it. String(take!(buf)) will be empty.
        # show(buf, Residue('A')) # This would write "A" to an empty buffer.

        show(buf, GAP) # GAP is Residue('-')
        @test String(take!(buf)) == "-"

        show(buf, XAA) # XAA is Residue('X')
        @test String(take!(buf)) == "X"
        
        # Test with invalid residues
        # The problem states: reinterpret(Residue, 0) -> "�"
        # Residue is a primitive type of 8 bits. So, UInt8(0) is a valid underlying value.
        # reinterpret(Residue, UInt8(0)) is how one might get such a Residue.
        # The character � (U+FFFD) is the Unicode replacement character.
        # This implies that Char(reinterpret(Residue, UInt8(0))) should yield something that print turns into "�".
        # This typically happens if the character is invalid for the IO context.
        # However, the current MIToS implementation of Char(res::Residue) returns 'X' for invalid residues.
        # If the test ` @test String(take!(buf)) == "�" ` is to pass,
        # then `Char(reinterpret(Residue, UInt8(0)))` must result in a character that serializes as "�".
        # For now, I will write the test as specified in the prompt.

        show(buf, reinterpret(Residue, UInt8(0))) # Invalid residue (0 is not in _RESIDUE_MAP keys)
        @test String(take!(buf)) == "�"
        
        show(buf, reinterpret(Residue, UInt8(25))) # Invalid residue (25 is not in _RESIDUE_MAP keys)
        @test String(take!(buf)) == "�"

        show(buf, reinterpret(Residue, UInt8(250))) # Another invalid residue
        @test String(take!(buf)) == "�"
    end

    @testset "Base.print(io::IO, res::Residue)" begin
        buf = IOBuffer()
        print(buf, Residue('A'))
        @test String(take!(buf)) == "A"

        print(buf, GAP)
        @test String(take!(buf)) == "-"

        print(buf, XAA)
        @test String(take!(buf)) == "X"

        # For print, the problem states the expected output for invalid is 'X'
        print(buf, reinterpret(Residue, UInt8(0))) # Invalid
        @test String(take!(buf)) == "X"

        print(buf, reinterpret(Residue, UInt8(25))) # Invalid
        @test String(take!(buf)) == "X"
        
        print(buf, reinterpret(Residue, UInt8(250))) # Invalid
        @test String(take!(buf)) == "X"
    end
end
