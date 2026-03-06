module GrassmannSymbolics

# ── source files ──────────────────────────────────────────────────────────────
include("generator.jl")
include("gexpr.jl")
include("berezin.jl")
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
export is_even
export is_odd
export is_grassmann

end # module GrassmannSymbolics
