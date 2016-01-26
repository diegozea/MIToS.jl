# using Base.Test
using IndexedArrays

const pdb_list = ["2HHB", "1IRD", "1DN3"]

print("""

Tests for IndexedArray
=======================
""")

print("""

Creation
--------
""")
const pdbs = IndexedArray(pdb_list)
@test pdbs.items == pdb_list
@test_throws IndexedArrayError IndexedArray([1,2,3,3])

print("""

Mapping
-------
""")

print("""
Test selectvalue
""")
for i in 1:3
  @test pdbs[i] == pdb_list[i]
end
@test_throws BoundsError pdbs[4]

@test pdbs[ [2,3] ]  == ["1IRD", "1DN3"]
@test pdbs[ Bool[false, true, true] ] == ["1IRD", "1DN3"]
@test pdbs[ [1,2,3] .> 1 ] == ["1IRD", "1DN3"]

@test pdbs[:] == pdb_list

print("""
Test selectindex
""")
i = 1
for pdb in pdb_list
  @test findfirst(pdbs, pdb) == i
  i += 1
end

@test_throws KeyError findfirst(pdbs, "2TRX")

@test [ findfirst(pdbs, i) for i in ["2HHB", "1IRD"] ] == Int[1,2]

@test [ findfirst(pdbs, i) for i in pdbs ] == Int[1,2,3]

print("""

Comparisons and getindex
------------------------
""")

@test IndexedArray(pdbs[[2,3]]) == IndexedArray(["1IRD", "1DN3"])

@test IndexedArray(pdbs[[false, true, true]]) == IndexedArray(["1IRD", "1DN3"])
@test IndexedArray(pdbs[[true, true, false]]) != IndexedArray(["1IRD", "1DN3"])

print("""

Iterations
----------
""")

@test collect(pdbs) == ["2HHB", "1IRD", "1DN3"]

i = 1
for element in pdbs
  @test element == pdb_list[i]
  i += 1
end

print("""

swap! and copy
--------------
""")

const copyofpdbs = copy(pdbs)

swap!(copyofpdbs, 2, 3)
@test findfirst(copyofpdbs, "1IRD") == 3
@test findfirst(copyofpdbs, "1DN3") == 2

@test copyofpdbs[2] == pdbs[3]
@test copyofpdbs[3] == pdbs[2]
@test findfirst(copyofpdbs, "1IRD") == findfirst(pdbs, "1DN3")
@test findfirst(copyofpdbs, "1DN3") == findfirst(pdbs, "1IRD")

print("""

push!
-----
""")
push!(pdbs, "2TRX")
@test_throws IndexedArrayError push!(pdbs, "2TRX")
@test pdbs[4] == "2TRX"
@test findfirst(pdbs, "2TRX") == 4
