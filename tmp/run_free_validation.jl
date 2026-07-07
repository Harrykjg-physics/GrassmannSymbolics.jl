using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
include(joinpath(@__DIR__, "test", "runtests.jl"))
include(joinpath(@__DIR__, "examples", "Free_Wilson_and_Staggered.jl"))
validate_free_tensors()
write(joinpath(@__DIR__, "success.sentinel"), "ok")