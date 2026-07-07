using Pkg

slash = Char(47)
dot = Char(46)
root = @__DIR__() * string(slash, dot, dot)
Pkg.activate(root)

example_path = joinpath(root, string(:examples), string(:NJL, dot, :jl))
include(example_path)

result = validate_njl_tensor()

println(string(:NJL_validation_passed))
println(string(:candidate_entries, :_, length(result.candidates)))

sentinel_path = joinpath(root, string(:success, dot, :sentinel))
open(sentinel_path, string(:w)) do io
    println(io, string(:ok))
end
