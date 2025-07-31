using BenchmarkTools
using NESSie

const SUITE = BenchmarkGroup()

include("bem/matvecmul.jl")
include("bem/solve.jl")
include("bem/post.jl")
