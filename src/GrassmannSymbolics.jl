module GrassmannSymbolics

import Symbolics

# ── source files ──────────────────────────────────────────────────────────────
include("generator.jl")
include("gexpr.jl")
include("berezin.jl")
include("localtensor.jl")
include("coefficientarray.jl")
include("show.jl")

# ── public API ────────────────────────────────────────────────────────────────
export GrassmannGenerator
export GrassmannExpr
export BerezinMeasure

export grassmann          # grassmann(:η), grassmann(:η, 1:3)
export dual               # dual(η) → η̄
export @grassmann         # @grassmann η θ ψ

export d                  # d(η) → BerezinMeasure

export scalar_part
export grassmann_part
export generators
export get_coeff
export map_coefficients
export map_generators
export substitute_generators
export simplify_coefficients
export substitute_coefficients
export symbolically_equal
export is_even
export is_odd
export is_grassmann

export integrate
export contract_grassmann
export local_grassmann_tensor
export LocalChannel
export NearestNeighborTensorSpec
export split_hopping_channel
export compile_nearest_neighbor_tensor
export LocalActionTerm
export OnsiteTerm
export FactorizedHoppingTerm
export NearestNeighborAction
export nearest_neighbor_tensor_spec
export coefficient_array
export coefficient_tensor
export coefficient_array_shape
export occupation_bit_tuples
export occupation_bits_from_index
export occupation_bits_from_indices
export occupation_index_from_bits
export occupied_legs_from_bits
export parameter_substitutions

end # module GrassmannSymbolics
