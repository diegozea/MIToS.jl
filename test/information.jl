using Base.Test
using MIToS.Information

## Fixed Pseudocount ##

# Constructor
none = zero(Fixed)
@test none.λ == zero(Float64)

# Copy & deepcopy
cnone = copy(none)
@test cnone.λ == zero(Float64)
dcnone = deepcopy(none)
@test dcnone.λ == zero(Float64)

## Pseudofrequencies ##

none = zeros(Pseudofrequencies)
@test none.α == zero(Float64)
@test none.β == zero(Float64)
@test none.Gab == zeros(Float64,(20,20))

# Copy & deepcopy
cnone = copy(none)
dcnone = deepcopy(none)

none.Gab[1,1] = 1.0

@test cnone.α == zero(Float64)
@test cnone.β == zero(Float64)
@test cnone.Gab[1,1] == zero(Float64)

@test dcnone.α == zero(Float64)
@test dcnone.β == zero(Float64)
@test dcnone.Gab[1,1] == zero(Float64)

## ResidueProbabilities ##

none = zeros(ResidueProbabilities)
@test none.Pa == zeros(Float64, 20)

# Copy & deepcopy
cnone = copy(none)
dcnone = deepcopy(none)

none.Pa[1] = 1.0

@test cnone.Pa[1] == zero(Float64)
@test dcnone.Pa[1] == zero(Float64)

## ResidueProbabilities ##

none =  zeros(ResiduePairProbabilities)
@test none.Pab == zeros(Float64, (20,20))
@test none.Pa == zeros(Float64, 20)
@test none.Pb == zeros(Float64, 20)

# Copy & deepcopy
cnone = copy(none)
dcnone = deepcopy(none)

none.Pab[1,1] = 1.0
none.Pa[1] = 1.0
none.Pb[1] = 1.0

@test cnone.Pab[1,1] == zero(Float64)
@test cnone.Pa[1] == zero(Float64)
@test cnone.Pb[1] == zero(Float64)

@test dcnone.Pab[1,1] == zero(Float64)
@test dcnone.Pa[1] == zero(Float64)
@test dcnone.Pb[1] == zero(Float64)

