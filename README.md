# GrassmannSymbolics.jl

[![CI](https://github.com/Harrykjg-physics/GrassmannSymbolics/actions/workflows/CI.yml/badge.svg)](https://github.com/Harrykjg-physics/GrassmannSymbolics/actions/workflows/CI.yml)
[![Documentation](https://github.com/Harrykjg-physics/GrassmannSymbolics/actions/workflows/Documentation.yml/badge.svg)](https://github.com/Harrykjg-physics/GrassmannSymbolics/actions/workflows/Documentation.yml)

`GrassmannSymbolics.jl` is a Julia package for symbolic calculations with
Grassmann, or anticommuting, variables.  It is designed for deriving local
Grassmann tensor-network tensors from fermionic lattice actions with
nearest-neighbor interactions.

The public documentation site is:

```text
https://harrykjg-physics.github.io/GrassmannSymbolics/
```

The source repository is:

```text
https://github.com/Harrykjg-physics/GrassmannSymbolics
```

## Features

- Grassmann generators with indexed and dual/barred variants.
- Sparse exterior-algebra expressions with automatic anticommutation signs.
- Nilpotent finite exponentials.
- Berezin integration with explicit measure order.
- Symbolic coefficients via `Symbolics.jl`.
- Local Grassmann tensor construction from onsite actions and factorized
  nearest-neighbor hopping channels.
- Coefficient tensor extraction as ordinary Julia `Array` objects.
- Model examples for Wilson/staggered fermions, Gross-Neveu-Wilson models,
  Schwinger model, two-color QCD, NJL, and 1D Hubbard.

## Installation

Before registration in Julia's General registry, install from the GitHub
repository:

```julia
using Pkg
Pkg.add(url="https://github.com/Harrykjg-physics/GrassmannSymbolics")
```

For local development:

```julia
using Pkg
Pkg.develop(path="/path/to/GrassmannSymbolics")
```

After registration in General:

```julia
using Pkg
Pkg.add("GrassmannSymbolics")
```

## Quick start

```julia
using GrassmannSymbolics

@grassmann eta xi

eta * xi      # eta*xi
xi * eta      # -eta*xi
eta * eta     # 0

exp(eta * xi) # 1 + eta*xi

integrate(eta * xi, [xi, eta]) # 1
integrate(eta * xi, [eta, xi]) # -1
```

With symbolic coefficients:

```julia
using Symbolics
@variables m

@grassmann psi
psibar = dual(psi)

expr = exp(m * psibar * psi)
```

## Local tensor kernel

```julia
@grassmann chi eta
chibar = dual(chi)
etabar = dual(eta)

T = local_grassmann_tensor(
    3chibar * chi,
    [-chibar * eta + (1 // 2) * etabar * chi],
    [chi, chibar],
)
```

The physical variables `chi, chibar` are integrated out, while the auxiliary
variables `eta, etabar` remain as tensor legs.

## Coefficient arrays

Grassmann tensors can be materialized as Julia arrays after choosing leg groups
and parameter substitutions:

```julia
@grassmann a b c d
expr = 2a*b + 3c*d + 5a*b*c*d

groups = [[a, b], [c, d]]
A = coefficient_array(expr, groups)
```

By default grouped indices use a parity-preserving order:

```text
(0, 0) -> 1
(1, 1) -> 2
(0, 1) -> 3
(1, 0) -> 4
```

## Examples

The `examples/` directory contains model derivations and validation helpers:

```text
Free_Wilson_and_Staggered.jl
Simple_quardratic_model.jl
Single_Flavor_Gross_Neveu_Wilson.jl
Schwinger_model_theta_term.jl
two_color_QCD.jl
two_flavor_Gross_Neveu_Wilson.jl
three_flavor_Gross_Neveu_Wilson.jl
NJL.jl
1D_Hubbard.jl
coefficient_tensors.jl
```

The detailed derivation workflow is documented in:

```text
docs/grassmann_tensor_derivation.tex
output/pdf/grassmann_tensor_derivation.pdf
```

## Two-flavor GNW erratum

In `Ref/2d_gn.pdf`, the impurity derivation for the two-flavor
Gross-Neveu-Wilson model has a sign swap in Eq. (3.79) and Eq. (3.80).  The
direct symbolic Berezin insertion gives:

```julia
b116 = (m1 + 2 - im*H) * (m2 + 2 + im*H) + (gσ2 - gπ2)/2
c116 = (m2 + 2 - im*H) * (m2 + 2 + im*H) + (gσ2 + gπ2)/2
```

The PDF writes the opposite `gπ2` signs.  Correcting these two coefficients
makes the pure tensor and the chiral, pseudoscalar-singlet, and
pseudoscalar-triplet impurity tensors agree with direct symbolic integration to
machine precision.

## Running tests

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

## Building documentation locally

```julia
using Pkg
Pkg.activate("docs")
Pkg.develop(path="..")
Pkg.instantiate()
include("docs/make.jl")
```

## Julia package registration

The package metadata is already present in `Project.toml`:

```toml
name = "GrassmannSymbolics"
uuid = "4c7e4f2a-8b3d-4e9a-b5c1-3f2d1e0a9c7b"
version = "0.1.0"
```

The standard registration route is JuliaRegistrator:

1. Push the release commit to GitHub.
2. On the release commit, issue, or pull request, comment:

   ```text
   @JuliaRegistrator register
   ```

3. Review the generated pull request in
   `https://github.com/JuliaRegistries/General`.
4. After the registry pull request is merged, TagBot can create the matching
   GitHub tag and release from the registered version.

See also `docs/src/registration.md`.

## License

This package is distributed under the MIT license.  See `LICENSE`.
