# using Base.Test
# using MIToS.Utils

print("""

Test General Utils
==================
""")

print("""
deleteitems!
""")
vec = collect("Hola mundo\n")
deleteitems!(vec, Char[' ', '\n'])
@test vec == Char['H', 'o', 'l', 'a', 'm', 'u', 'n', 'd', 'o']

vec = Int[ char for char in "Hola mundo\n" ]
deleteitems!(vec, IntSet(Int[' ', '\n']) )
@test vec == Int['H', 'o', 'l', 'a', 'm', 'u', 'n', 'd', 'o']

print("""
get_n_words
""")
line = "#=GF AC PF00571"
line_n = "#=GF AC PF00571\n"
@test get_n_words(line, 1) == ASCIIString[line]
@test get_n_words(line, 2) == ASCIIString["#=GF", "AC PF00571"]
@test get_n_words(line, 3) == ASCIIString["#=GF", "AC", "PF00571"]
@test get_n_words(line, 4) == ASCIIString["#=GF", "AC", "PF00571"]
for i in 1:4
  @test get_n_words(line, i) == get_n_words(line_n, i)
end

@test get_n_words("\n",1) == ASCIIString[]
@test get_n_words("#", 1) == ASCIIString["#"]

let str = "#=GF CC   (Römling U.  and Galperin M.Y. “Bacterial cellulose\n"

  @test get_n_words(str, 3) == UTF8String["#=GF", "CC", "(Römling U.  and Galperin M.Y. “Bacterial cellulose"]
end

let str = "#=GF CC   not present in all SecA2–SecY2 systems. This family of Asp5 is\n"

  @test get_n_words(str, 3) == UTF8String["#=GF", "CC", "not present in all SecA2–SecY2 systems. This family of Asp5 is"]
end

print("""
hascoordinates
""")
@test hascoordinates("O83071/192-246")
@test !hascoordinates("O83071")

print("""
select_element
""")
@test select_element([1]) == 1
# print("Test warn: ")
# @test select_element([1,2]) == 1
@test_throws ErrorException select_element([])

print("""

matrix2list
-----------
""")

let mat = [ 1 2 3
            4 5 6
            7 8 9 ]

  @test matrix2list(mat) == [2, 3, 6]
  @test matrix2list(mat, diagonal=true) == [1, 2, 3, 5, 6, 9]
  @test matrix2list(mat, part="lower") == [4, 7, 8]
  @test matrix2list(mat, part="lower", diagonal=true) == [1, 4, 7, 5, 8, 9]
end

print("""

list2matrix
-----------
""")

let mat = [ 1 2 3
            2 5 6
            3 6 9 ]

  @test triu(list2matrix([2, 3, 6], 3), 1) == triu(mat, 1)
  @test list2matrix([1, 2, 3, 5, 6, 9], 3, diagonal=true) == mat
end

print("""

Test findobjects
================
""")

import Base: ==, hash

using AutoHashEquals

@auto_hash_equals immutable Dummy
  string::ASCIIString
  int::Int
end

print("""
findobjects, isobject and AbstractTest
""")

@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Not(Is(:int,4))) == [1, 2]
@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Is(:int,4)) == [3]

@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Is(:int, x -> x > 2 )) == [2, 3]
@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Not(Is(:int, x -> x > 2 ))) == [1]

@test findobjects([ Dummy("abc", 2), Dummy("abcd", 2), Dummy("bc", 2) ], Is(:string, r"^ab")) == [1, 2]
@test findobjects([ Dummy("abc", 2), Dummy("abcd", 2), Dummy("bc", 2) ], Not(Is(:string, r"^ab"))) == [3]
@test findobjects([ Dummy("abc", 2), Dummy("abcd", 2), Dummy("bc", 2) ], Is(:string, r"^nothere$")) == []

@test findobjects([ Dummy("H", 2), Dummy("C", 2), Dummy("O", 2) ], Not(Is(:string, "H"))) == [2, 3]

@test findobjects([ Dummy("H", 2), Dummy("C", 2), Dummy("O", 2) ], In(:string, ["C", "O"])) == [2, 3]
@test findobjects([ Dummy("H", 2), Dummy("C", 2), Dummy("O", 2) ], Not(In(:string, ["C", "O"]))) == [1]

print("""
collectobjects & isobject
""")

@test collectobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Is(:int,4)) == [Dummy("c", 4)]
@test collectobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Is(:int, x -> x > 2 )) == [Dummy("b", 3), Dummy("c", 4)]
@test collectobjects([ Dummy("abc", 2), Dummy("abcd", 2), Dummy("bc", 2) ], Not(Is(:string, r"^ab"))) == [Dummy("bc", 2)]
@test collectobjects([ Dummy("abc", 2), Dummy("abcd", 2), Dummy("bc", 2) ], Is(:string, r"^nothere$")) == []

print("""
collectcaptures & capture
""")

@test map(isnull, collectcaptures([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], :string, Is( :int, x -> x > 2 )))  == [true, false, false]
@test map(isnull, collectcaptures([ "NotDummy", 42, Dummy("c", 4) ], :string, Is( :int, x -> x > 2 )))   == [true, true, false]

print("""
collectcaptures & guess_type
""")

@test isa(collectcaptures([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], :string, Is( :int, 3 )), Array{Nullable{ASCIIString},1})
@test isa(collectcaptures([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], :int, Is( :int, 3 )), Array{Nullable{Int},1})
