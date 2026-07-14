# Coefficient arrays

Grassmann tensor expressions are sparse symbolic objects.  Tensor-network codes
often need dense or sparse ordinary Julia arrays after choosing:

1. numerical parameter values;
2. tensor-leg grouping;
3. the occupation-bit order inside each grouped index.

Use:

```julia
coefficient_array(tensor, leg_groups; substitutions, index_order=:parity)
```

or the alias:

```julia
coefficient_tensor(...)
```

## Parity-preserving grouped indices

The default `index_order=:parity` puts even states first and odd states second.
For two bits:

```text
(0, 0) -> 1
(1, 1) -> 2
(0, 1) -> 3
(1, 0) -> 4
```

For four bits, the same rule gives an index dimension of length `16`, with all
even bitstrings before odd bitstrings.

## Example

```julia
@grassmann a b c d
expr = 2a*b + 3c*d + 5a*b*c*d
groups = [[a, b], [c, d]]

A = coefficient_array(expr, groups; element_type=Int)
size(A) == (4, 4)
```

## Parameter substitutions

For models with a `params` or `parameters` named tuple:

```julia
rules = parameter_substitutions(model.params, (; M1=0.1, M2=0.2))
array = coefficient_array(model.tensor, groups; substitutions=rules)
```

The helper functions in `examples/coefficient_tensors.jl` wrap this pattern for
the example models.
