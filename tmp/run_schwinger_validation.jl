import Pkg

Pkg.activate(@__DIR__)
include(joinpath(
    @__DIR__,
    string(:examples),
    string(:Schwinger_model_theta_term, Char(46), :jl),
))
validate_schwinger_fermion_tensor()

open(
    joinpath(@__DIR__, string(:success, Char(46), :sentinel));
    write=true,
    create=true,
    truncate=true,
) do io
    write(io, string(:ok))
end
