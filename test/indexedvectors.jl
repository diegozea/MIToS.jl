using Base.Test
using MIToS.MSA

const pdb_list = ["2HHB", "1IRD", "1DN3"]

print("""

Tests for IndexedVector
=======================
""")

print("""

Creation
--------
""")
const pdbs = IndexedVector(pdb_list)
@test pdbs.values == pdb_list
# @test_throws "3 is more than one time in the vector." IndexedVector([1,2,3,3])

print("""

Mapping
-------
""")

print("""
Test selectvalue
""")
for i in 1:3
  @test selectvalue(pdbs, i) == pdb_list[i]
end
@test_throws BoundsError selectvalue(pdbs, 4)

@test selectvalue(pdbs, [2,3]) == ["1IRD", "1DN3"]
@test selectvalue(pdbs, Bool[false, true, true]) == ["1IRD", "1DN3"]
@test selectvalue(pdbs, [1,2,3] .> 1 ) == ["1IRD", "1DN3"]

@test selectvalue(pdbs) == pdb_list

print("""
Test selectindex
""")
i = 1
for pdb in pdb_list
  @test selectindex(pdbs, pdb) == i
  i += 1
end

@test_throws KeyError selectindex(pdbs, "2TRX")

@test selectindex(pdbs, ["2HHB", "1IRD"]) == Int[1,2]

@test selectindex(pdbs) == Int[1,2,3]

print("""

Comparisons and getindex
------------------------
""")

@test pdbs[[2,3]] == IndexedVector(["1IRD", "1DN3"])

@test pdbs[[false, true, true]] == IndexedVector(["1IRD", "1DN3"])
@test pdbs[[true, true, false]] != IndexedVector(["1IRD", "1DN3"])

print("""

Iterations
----------
""")

@test collect(pdbs) == Tuple{ASCIIString,Int}[("2HHB",1), ("1IRD",2), ("1DN3",3)]

i = 1
for element in pdbs
  @test isa(element, Tuple{ASCIIString,Int})
  @test element[1] == pdb_list[i]
  @test element[2] == i
  i += 1
end

print("""

swap! and copy
--------------
""")

const copyofpdbs = copy(pdbs)

swap!(copyofpdbs, 2, 3)
@test selectindex(copyofpdbs, "1IRD") == 3
@test selectindex(copyofpdbs, "1DN3") == 2

@test selectvalue(copyofpdbs, 2) == selectvalue(pdbs, 3)
@test selectvalue(copyofpdbs, 3) == selectvalue(pdbs, 2)
@test selectindex(copyofpdbs, "1IRD") == selectindex(pdbs, "1DN3")
@test selectindex(copyofpdbs, "1DN3") == selectindex(pdbs, "1IRD")

print("""

push!
-----
""")
push!(pdbs, "2TRX")
@test selectvalue(pdbs, 4) == "2TRX"
@test selectindex(pdbs, "2TRX") == 4




