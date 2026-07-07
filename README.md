# GrassmannSymbolics.jl

A Julia package for symbolic computation with **Grassmann (anticommuting) numbers** and the exterior algebra.  
Inspired by the Python package [grassmanntn](https://github.com/ayosprakob/grassmanntn) and leveraging [Symbolics.jl](https://symbolics.juliasymbolics.org/) for symbolic coefficients.

## Features

- **Grassmann generators** with optional subscript indices (η₁, η₂, …) and conjugates (η̄)
- **Full exterior algebra arithmetic**: addition, subtraction, scalar multiplication/division, anticommutative products
- **Automatic nilpotency**: η² = 0, higher products terminate automatically
- **`exp` function**: finite Taylor series via nilpotency
- **Berezin integration**: composable `d(η)` operators
- **Symbolic coefficients**: works seamlessly with `Symbolics.jl` variables
- **Unicode display**: pretty-prints with subscripts and combining overbar

## Installation

```julia
# From the Julia REPL (package mode):
pkg> dev path/to/GrassmannSymbolics
# or, once registered:
pkg> add GrassmannSymbolics
```

## Quick Start

```julia
using GrassmannSymbolics

# Create generators
@grassmann η θ ψ          # η, θ, ψ
ηbar = dual(η)             # η̄
η1, η2 = grassmann(:η, 1:2)  # η₁, η₂

# Anticommutation
η * θ                      # ηθ
θ * η                      # -ηθ
η * η                      # 0   (nilpotency)

# Arithmetic
2η + 3θ                    # 2η + 3θ
(1 + η) * (1 + θ)          # 1 + θ + η + ηθ

# Powers and exp
η^2                        # 0
exp(η * θ)                 # 1 + ηθ
exp(1 + η)                 # ℯ + ℯ·η
exp(η + θ)                 # 1 + η + θ  because (η + θ)² = 0

# Berezin integration
d(η) * GrassmannExpr(η)    # 1
d(η) * (η * θ)             # θ
d(θ) * d(η) * (η * θ)     # 1
d(η) * d(θ) * (η * θ)     # -1

# With Symbolics.jl
using Symbolics
@variables a b m
a*η + b*θ                  # symbolic coefficients
exp(m * η * θ)             # 1 + m·ηθ
```

## API Reference

| Name | Description |
|------|-------------|
| `@grassmann η θ ...` | Batch-create generators bound to local variables |
| `grassmann(:η)` | Create a single generator |
| `grassmann(:η, 1:n)` | Create a vector of indexed generators η₁…ηₙ |
| `dual(η)` | Return the conjugate η̄ |
| `d(η)` | Create a Berezin measure; chain with `*` |
| `scalar_part(e)` | Extract the grade-0 coefficient |
| `grassmann_part(e)` | Extract all grade ≥ 1 terms |
| `get_coeff(e, [η,θ])` | Extract the coefficient of a given monomial |
| `map_coefficients(f, e)` | Transform scalar coefficients without changing monomials |
| `map_generators(f, e)` | Relabel tensor legs and restore canonical signs |
| `substitute_generators(e, rules)` | Instantiate local tensor legs on lattice sites |
| `simplify_coefficients(e)` | Simplify all coefficients with Symbolics.jl |
| `substitute_coefficients(e, rules)` | Substitute Symbolics.jl variables coefficient-wise |
| `symbolically_equal(a, b)` | Compare expressions after coefficient simplification |
| `integrate(e, order)` | Perform a Berezin integral with an explicit measure order |
| `contract_grassmann(tensors, pairs)` | Contract oriented legs with Gaussian link measures |
| `local_grassmann_tensor(...)` | Integrate physical site fields and leave auxiliary tensor legs |
| `generators(e)` | List all generators appearing in `e` |
| `is_even(e)` | True if all terms have even grade |
| `is_odd(e)` | True if all terms have odd grade |
| `is_grassmann(e)` | True if any term has grade ≥ 1 |

## Local Grassmann tensors

The model-independent site kernel computes

```math
\mathcal T = \int D\psi\,e^{-S_{\mathrm{onsite}}}\prod_k e^{B_k},
```

where `B_k` are the local auxiliary-field factors obtained from directed
nearest-neighbor hopping terms:

```julia
@grassmann ψ η
ψbar, ηbar = dual(ψ), dual(η)

T = local_grassmann_tensor(
    3ψbar * ψ,
    [-ψbar * η + (1 // 2) * ηbar * ψ],
    [ψ, ψbar],
)
# T = -3 - 1/2 ηηbar
```

The derivation conventions and the model-by-model audit are documented in
[`docs/grassmann_tensor_derivation.tex`](docs/grassmann_tensor_derivation.tex).
## Running Tests

```julia
using Pkg
Pkg.test("GrassmannSymbolics")
```

## Design Notes

- **Representation**: `GrassmannExpr` stores a sparse `Dict{Vector{GrassmannGenerator}, Any}` mapping each sorted monomial to its coefficient.  The empty vector `[]` represents the scalar "1".
- **Sign computation**: arbitrary input lists are normalized by inversion counting; products of already-canonical monomials use a linear-time ordered merge that tracks permutation parity.
- **exp**: uses the decomposition `exp(c + X) = exp(c) · Σ Xᵖ/p!`, which terminates because `X` is nilpotent.
- **Berezin integral**: applies `∫dg` to each term by locating `g` in the sorted monomial, extracting it with the appropriate sign `(-1)^(position-1)`.
