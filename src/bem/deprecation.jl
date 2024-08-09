# removed in v1.4
@deprecate solve_implicit(lt, model) solve(lt, model; method=:gmres)
