using ParallelTestRunner
import NESSie

const init_code = quote
    using NESSie

    include(nessie_data_path("../test/testsetup.jl"))
end

testsuite = find_tests(@__DIR__)
delete!(testsuite, "testsetup")

# Compat: remove file extension from test names
replace!(e -> endswith(e, ".jl") ? e[1:end-3] : e, ARGS)

runtests(NESSie, ["--verbose", ARGS...]; init_code, testsuite)
