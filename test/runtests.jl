using Test

include("tests.jl")

# Doctests
if VERSION >= v"1.7.0" # Julia 1.7 changed the default random number generator Mersenne Twister to Xoshiro256++
    doctest(MIToS)
end

print("""

----- =D -----

""")
