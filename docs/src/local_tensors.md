# Local tensor compiler

The central tensor kernel is:

```julia
local_grassmann_tensor(onsite_action, bond_exponents, measure_order)
```

It computes

```math
\mathcal T =
\int D\psi\,
\exp(-S_{\mathrm{onsite}})
\prod_k \exp(B_k),
```

where:

- `onsite_action` is the local action involving physical fields;
- `bond_exponents` are local auxiliary-field factors from hopping terms;
- `measure_order` fixes the Berezin convention for integrating the physical
  fields.

## Structured nearest-neighbor actions

For factorized nearest-neighbor terms, use:

```julia
@grassmann χ η
χbar = dual(χ)
ηbar = dual(η)

action = NearestNeighborAction(
    [χ, χbar];
    terms=[
        OnsiteTerm(3χbar * χ; tag=:mass),
        FactorizedHoppingTerm(
            χbar,
            η,
            χ;
            left_scale=-1,
            right_scale=1 // 2,
            tag=:hopping,
        ),
    ],
    leg_order=[η, ηbar],
)

compiled = compile_nearest_neighbor_tensor(action)
```

`compiled` contains:

- `tensor`: the local Grassmann tensor;
- `legs`: tensor-leg order;
- `channels`: auxiliary-channel metadata;
- `measure_order`: physical-field Berezin order;
- `metadata`: model-specific annotations.

Dense hopping matrices, gauge-field sums, and four-fermion decompositions still
need model-specific preprocessing before they become factorized hopping terms.
