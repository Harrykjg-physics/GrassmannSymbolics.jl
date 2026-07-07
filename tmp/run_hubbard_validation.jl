using Pkg

slash = Char(47)
dot = Char(46)
root = @__DIR__() * string(slash, dot, dot)
Pkg.activate(root)

example_path = joinpath(root, string(:examples), string(1, :D_Hubbard, dot, :jl))
include(example_path)

model = validate_hubbard_tensor()

println(string(:Hubbard_validation_passed))
println(string(:tensor_terms, :_, length(model.tensor.terms)))

sentinel_path = joinpath(root, string(:success, dot, :sentinel))
open(sentinel_path, string(:w)) do io
    println(io, string(:ok))
end
