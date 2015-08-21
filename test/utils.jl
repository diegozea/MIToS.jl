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

print("""
hascoordinates
""")
@test hascoordinates("O83071/192-246")
@test !hascoordinates("O83071")

print("""
select_element
""")
@test select_element([1]) == 1
@test select_element([1,2]) == 1
@test_throws ErrorException select_element([])
