# Model examples

The `examples/` directory contains model-specific derivations and validation
scripts.  These examples show how the model-independent Grassmann kernel can be
used to reproduce published tensor-network formulas.

## Available models

```text
Free_Wilson_and_Staggered.jl
Wilson_Majorana_Fermion.jl
Two_Flavor_Staggered_Gross_Neveu.jl
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

## Wilson-Majorana and two-flavor staggered GN

`Wilson_Majorana_Fermion.jl` reconstructs Eq. (2.9) of
`Ref/Wilson_Majorana_Fermion_and_Two_Flavor_Staggered_GN.pdf` by direct
Berezin integration over the two Majorana spinors.  The four tensor legs are
the 4-bit super-indices `(i, j, k, l)`, so the coefficient array has shape
`(16, 16, 16, 16)`.

`Two_Flavor_Staggered_Gross_Neveu.jl` reconstructs Eq. (5.2) of the same
reference by integrating the four staggered components `psi_1,...,psi_4` with
the onsite interaction `exp(U psi_1 psi_2 psi_3 psi_4)`.  It uses the same
four 4-bit leg groups and therefore the same coefficient-array shape.

## Coefficient-tensor helpers

`examples/coefficient_tensors.jl` provides convenience constructors for turning
the symbolic tensors into Julia `Array` objects.  Large models can be queried
with `materialize=false` to inspect the shape and leg grouping without building
the full dense array.

For the two-flavor Gross-Neveu-Wilson model, the physical-parameter helper maps:

```julia
(M1, M2, gs, gpi, r1, r2, mu1, mu2, H)
```

to the symbolic parameters used in the derivation:

```julia
M1, M2, gs2, gp2, H, q, z1, z2
```

where `q = 1/sqrt(2)` and `z_f = exp(mu_f/2)` by default.

## Two-flavor GNW erratum

The reference `Ref/2d_gn.pdf` contains a sign swap in the impurity derivation for
the two-flavor Gross-Neveu-Wilson model.  In Eq. (3.79) and Eq. (3.80), the
`g_pi^2` signs in the coefficients corresponding to `b116` and `c116` are
interchanged.

The direct symbolic Berezin insertion gives:

```julia
b116 = (m1 + 2 - im*H) * (m2 + 2 + im*H) + (gσ2 - gπ2)/2
c116 = (m2 + 2 - im*H) * (m2 + 2 + im*H) + (gσ2 + gπ2)/2
```

The PDF writes the opposite signs.  Numerically, using the corrected pair makes
the pure tensor and the chiral, pseudoscalar-singlet, and pseudoscalar-triplet
impurity tensors agree with direct symbolic integration to machine precision.
