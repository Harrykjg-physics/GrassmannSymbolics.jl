using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
include(joinpath(@__DIR__, "tmp", "diagnose_free.jl"))
write(joinpath(@__DIR__, "success.sentinel"), "ok")