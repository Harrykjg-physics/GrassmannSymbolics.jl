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

## Installation and package registration

This repository is a standard Julia package named `GrassmannSymbolics`.
For local development, register it in your active Julia environment with:

```julia
using Pkg
Pkg.develop(path="path/to/GrassmannSymbolics")
```

For read-only use directly from a local checkout, install it by path:

```julia
using Pkg
Pkg.add(path="path/to/GrassmannSymbolics")
```

When working inside this repository, activate the package project directly:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using GrassmannSymbolics
```

To publish the package to Julia's General registry later, push this repository
to a public Git host, tag a release matching `version` in `Project.toml`, and
open a registration request with JuliaRegistrator. The package metadata needed
for that workflow lives in `Project.toml`.

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
| `LocalChannel(exponent; legs, tag)` | One local factor `exp(exponent)` plus the tensor legs it introduces |
| `NearestNeighborTensorSpec(...)` | A reusable nearest-neighbor local tensor specification |
| `split_hopping_channel(left, aux, right; ...)` | Split a rank-one hopping term into two auxiliary Grassmann factors |
| `compile_nearest_neighbor_tensor(spec)` | Compile a nearest-neighbor tensor spec to `(tensor, legs, channels, ...)` |
| `OnsiteTerm(action; tag)` | Structured onsite contribution to a local action |
| `FactorizedHoppingTerm(left, aux, right; ...)` | Structured rank-one nearest-neighbor hopping contribution |
| `NearestNeighborAction(measure_order; terms, ...)` | Build a tensor spec from structured local action terms |
| `nearest_neighbor_tensor_spec(action)` | Lower a structured action to `NearestNeighborTensorSpec` |
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

For already-factorized nearest-neighbor hopping terms, the higher-level compiler
keeps the tensor, leg order, channel metadata, and measure convention bundled
together:

```julia
@grassmann chi eta
chibar, etabar = dual(chi), dual(eta)

channels = split_hopping_channel(
    chibar,
    eta,
    chi;
    left_scale=-1,
    right_scale=1 // 2,
    tag=:forward_hopping,
)

spec = NearestNeighborTensorSpec(
    3chibar * chi,
    [chi, chibar],
    channels;
    leg_order=[eta, etabar],
    metadata=Dict(:model => :one_component_test),
)

compiled = compile_nearest_neighbor_tensor(spec)
compiled.tensor
```

The same example can be written one level higher as a structured local action:

```julia
action = NearestNeighborAction(
    [chi, chibar];
    terms=[
        OnsiteTerm(3chibar * chi; tag=:mass),
        FactorizedHoppingTerm(
            chibar,
            eta,
            chi;
            left_scale=-1,
            right_scale=1 // 2,
            tag=:forward_hopping,
        ),
    ],
    leg_order=[eta, etabar],
)

compiled = compile_nearest_neighbor_tensor(action)
```

This compiler layer assumes each nearest-neighbor term has already been written
as local left/right factors. Dense hopping matrices, gauge-field sums, or
four-fermion channels still need a model-specific decomposition step before
calling `compile_nearest_neighbor_tensor`.

The derivation conventions and the model-by-model audit are documented in
[`docs/grassmann_tensor_derivation.tex`](docs/grassmann_tensor_derivation.tex).

## Running Tests

From this repository:

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

After `Pkg.develop(path="path/to/GrassmannSymbolics")` from another environment:

```julia
using Pkg
Pkg.test("GrassmannSymbolics")
```

## Design Notes

- **Representation**: `GrassmannExpr` stores a sparse `Dict{Vector{GrassmannGenerator}, Any}` mapping each sorted monomial to its coefficient.  The empty vector `[]` represents the scalar "1".
- **Sign computation**: arbitrary input lists are normalized by inversion counting; products of already-canonical monomials use a linear-time ordered merge that tracks permutation parity.
- **exp**: uses the decomposition `exp(c + X) = exp(c) · Σ Xᵖ/p!`, which terminates because `X` is nilpotent.
- **Berezin integral**: applies `∫dg` to each term by locating `g` in the sorted monomial, extracting it with the appropriate sign `(-1)^(position-1)`.
