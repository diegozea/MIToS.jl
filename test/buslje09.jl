using Base.Test
using MIToS.Information
using MIToS.Utils

const Gaoetal2011 = joinpath(pwd(), "data", "Gaoetal2011.fasta")

gao11_buslje09(measure) = joinpath(pwd(), "data", string("data_Gaoetal2011_soft_Busljeetal2009_measure_", measure, ".txt"))

## Column numbers for the output of Buslje et. al. 2009
const SCORE = 9
const ZSCORE = 12

let data = readdlm(gao11_buslje09("MI")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=false)
  @test_approx_eq_eps maximum(abs(convert(Vector{Float64}, data[:, SCORE]) .- matrix2list(results[2]))) 0.0 0.000001
  @test cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])) > 0.90

  println(cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])))
end

let data = readdlm(gao11_buslje09("MI_clustering")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=true, apc=false)
  @test_approx_eq_eps maximum(abs(convert(Vector{Float64}, data[:, SCORE]) .- matrix2list(results[2]))) 0.0 0.000001
  @test cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])) > 0.90

  println(cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])))
end

let data = readdlm(gao11_buslje09("MI_APC")); results = buslje09(Gaoetal2011, FASTA, lambda=0.0, clustering=false, apc=true)
  #@test_approx_eq_eps maximum(abs(convert(Vector{Float64}, data[:, SCORE]) .- matrix2list(results[2]))) 0.0 0.000001
  @test cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])) > 0.90

  println(cor(convert(Vector{Float64}, data[:, ZSCORE]), matrix2list(results[3])))
end

