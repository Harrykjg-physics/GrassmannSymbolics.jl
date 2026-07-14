# GrassmannSymbolics.jl

`GrassmannSymbolics.jl` is a Julia package for symbolic calculations in
Grassmann exterior algebras.  It was built to derive local Grassmann tensor
network tensors from fermionic lattice actions with nearest-neighbor hopping,
while keeping Berezin measure order, tensor-leg order, and coefficient-array
conventions explicit.

The package is useful when you need to:

- manipulate anticommuting generators symbolically;
- expand finite Grassmann exponentials;
- perform Berezin integration with a fixed measure convention;
- construct local tensor-network kernels from onsite actions and factorized
  nearest-neighbor hopping channels;
- materialize sparse Grassmann expressions as ordinary Julia `Array` coefficient
  tensors after substituting numerical model parameters;
- audit published Grassmann tensor formulas against direct symbolic
  computation.

## Documentation URL

After the documentation workflow is enabled on GitHub Pages, the public manual
is served at:

```text
https://harrykjg-physics.github.io/GrassmannSymbolics.jl/
```

## Installation

During development, install the package from a checkout:

```julia
using Pkg
Pkg.develop(path="/path/to/GrassmannSymbolics.jl")
```

After registration in Julia's General registry, users will be able to install it
with:

```julia
using Pkg
Pkg.add("GrassmannSymbolics")
```

## Package layout

```text
src/        Core Grassmann algebra, Berezin integration, tensor compiler.
examples/   Model derivations and coefficient-tensor helpers.
test/       Unit and model-validation tests.
docs/       This Documenter.jl manual and the longer LaTeX derivation note.
Ref/        Reference PDFs used for model-by-model audits.
```

## Main conventions

Grassmann monomials are stored in a canonical sorted order.  Whenever a user
requests a coefficient with a different generator order, the package applies the
corresponding anticommutation sign.  Berezin integration always takes an
explicit measure order; this avoids hiding sign conventions in global state.

Grouped coefficient-array indices use a parity-preserving order by default:

```text
(0, 0) -> 1
(1, 1) -> 2
(0, 1) -> 3
(1, 0) -> 4
```

The same idea extends to wider grouped indices: even parity states come first,
then odd parity states, each block in lexicographic order.

## Related derivation note

The model-independent derivation workflow is also available as a LaTeX note:

```text
docs/grassmann_tensor_derivation.tex
```

The compiled PDF artifact in this repository is:

```text
output/pdf/grassmann_tensor_derivation.pdf
```
