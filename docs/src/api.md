# API index

This page lists the main exported names.  See the manual pages for the workflow
context.

## Core types

- `GrassmannGenerator`
- `GrassmannExpr`
- `BerezinMeasure`

## Generator constructors

- `grassmann`
- `dual`
- `@grassmann`

## Algebra and predicates

- `scalar_part`
- `grassmann_part`
- `generators`
- `get_coeff`
- `map_coefficients`
- `map_generators`
- `substitute_generators`
- `simplify_coefficients`
- `substitute_coefficients`
- `symbolically_equal`
- `is_even`
- `is_odd`
- `is_grassmann`

## Berezin integration and tensor construction

- `d`
- `integrate`
- `contract_grassmann`
- `local_grassmann_tensor`
- `LocalChannel`
- `NearestNeighborTensorSpec`
- `split_hopping_channel`
- `compile_nearest_neighbor_tensor`
- `LocalActionTerm`
- `OnsiteTerm`
- `FactorizedHoppingTerm`
- `NearestNeighborAction`
- `nearest_neighbor_tensor_spec`

## Coefficient arrays

- `coefficient_array`
- `coefficient_tensor`
- `coefficient_array_shape`
- `occupation_bit_tuples`
- `occupation_bits_from_index`
- `occupation_bits_from_indices`
- `occupation_index_from_bits`
- `occupied_legs_from_bits`
- `parameter_substitutions`
