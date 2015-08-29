using Base.Test
using MIToS.Utils

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

Test findobjects
================
""")

immutable Dummy
  string::ASCIIString
  int::Int
end

@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Not(:int,4)) == [1, 2]
@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Is( :int,4)) == [3]

@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Is( :int, x -> x > 2 ))  == [2, 3]
@test findobjects([ Dummy("a", 2), Dummy("b", 3), Dummy("c", 4) ], Not(:int, x -> x > 2 ))  == [1]

@test findobjects([ Dummy("abc", 2), Dummy("abcd", 2), Dummy("bc", 2) ], Is( :string, r"^ab")) == [1, 2]
@test findobjects([ Dummy("abc", 2), Dummy("abcd", 2), Dummy("bc", 2) ], Not(:string, r"^ab")) == [3]

@test findobjects([ Dummy("H", 2), Dummy("C", 2), Dummy("O", 2) ], Not(:string, "H")) == [2, 3]
