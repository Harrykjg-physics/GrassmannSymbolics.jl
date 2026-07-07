using Pkg

slash = Char(47)
dot = Char(46)
root = @__DIR__() * string(slash, dot, dot)
Pkg.activate(root)

example_path = joinpath(root, string(:examples), string(:two_color_QCD, dot, :jl))
include(example_path)

model = validate_two_color_qcd_fermion_tensor()
source_model = derive_two_color_qcd_fermion_tensor(; diquark_source=true)

println(string(:two_color_QCD_validation_passed))
println(string(:fermion_tensor_terms, :_, length(model.tensor.terms)))
println(string(:diquark_source_tensor_terms, :_, length(source_model.tensor.terms)))

sentinel_path = joinpath(root, string(:success, dot, :sentinel))
open(sentinel_path, string(:w)) do io
    println(io, string(:ok))
end
