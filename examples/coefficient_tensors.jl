using GrassmannSymbolics
using Symbolics

include(joinpath(@__DIR__, "Free_Wilson_and_Staggered.jl"))
include(joinpath(@__DIR__, "Simple_quardratic_model.jl"))
include(joinpath(@__DIR__, "Single_Flavor_Gross_Neveu_Wilson.jl"))
include(joinpath(@__DIR__, "Schwinger_model_theta_term.jl"))
include(joinpath(@__DIR__, "two_color_QCD.jl"))
include(joinpath(@__DIR__, "two_flavor_Gross_Neveu_Wilson.jl"))
include(joinpath(@__DIR__, "three_flavor_Gross_Neveu_Wilson.jl"))
include(joinpath(@__DIR__, "NJL.jl"))
include(joinpath(@__DIR__, "1D_Hubbard.jl"))

# Unified coefficient Array constructors for the example models.
#
# Each function accepts either direct symbolic substitutions, e.g.
#   substitutions=Dict(model.parameters.m => 1.0)
# or a NamedTuple of parameter values when the model exposes `parameters` or
# `params`, e.g.
#   parameters=(; m=1.0, q=0.5)
#
# `materialize=false` returns a lightweight specification containing the source
# tensor/oracle, leg groups, and Array shape.  This is useful for large tensors.

function _model_parameters(model)
    hasproperty(model, :parameters) && return getproperty(model, :parameters)
    hasproperty(model, :params) && return getproperty(model, :params)
    return nothing
end

function _substitution_rules(model, parameter_values, substitutions)
    substitutions !== nothing && return substitutions
    params = _model_parameters(model)
    params === nothing && return nothing
    isempty(keys(parameter_values)) && return nothing
    return parameter_substitutions(params, parameter_values)
end

function _coefficient_array_result(source, groups; materialize::Bool=true,
        substitutions=nothing, simplify::Bool=false, element_type=Any,
        index_order::Symbol=:parity)
    if materialize
        return coefficient_array(source, groups;
            substitutions=substitutions,
            simplify=simplify,
            element_type=element_type,
            index_order=index_order,
        )
    end
    return (;
        source,
        leg_groups=groups,
        shape=coefficient_array_shape(groups),
        substitutions,
        simplify,
        element_type,
        index_order,
    )
end

_pair_groups(legs) = [collect(legs[1:2]), collect(legs[3:4]), collect(legs[5:6]), collect(legs[7:8])]
_four_groups_16(legs) = [collect(legs[1:4]), collect(legs[5:8]), collect(legs[9:12]), collect(legs[13:16])]

function free_wilson_coefficient_tensor(; parameters=(;), substitutions=nothing,
        published::Bool=false, materialize::Bool=true, simplify::Bool=false, element_type=Any)
    model = derive_wilson_tensor()
    tensor = published ? model.published_tensor : model.tensor
    groups = [collect(group) for group in values(model.legs)]
    rules = _substitution_rules(model, parameters, substitutions)
    return _coefficient_array_result(tensor, groups;
        materialize, substitutions=rules, simplify, element_type)
end

function free_staggered_coefficient_tensor(; substitutions=nothing,
        published::Bool=false, materialize::Bool=true, simplify::Bool=false, element_type=Any)
    model = derive_staggered_tensor()
    tensor = published ? model.published_tensor : model.tensor
    groups = [collect(group) for group in values(model.legs)]
    return _coefficient_array_result(tensor, groups;
        materialize, substitutions, simplify, element_type)
end

function simple_quadratic_coefficient_tensor(; parameters=(;), substitutions=nothing,
        materialize::Bool=true, simplify::Bool=false, element_type=Any)
    model = derive_simple_quadratic_tensor()
    rules = _substitution_rules(model, parameters, substitutions)
    return _coefficient_array_result(model.tensor, _pair_groups(model.legs);
        materialize, substitutions=rules, simplify, element_type)
end

function single_flavor_gnw_coefficient_tensor(; parameters=(;), substitutions=nothing,
        materialize::Bool=true, simplify::Bool=false, element_type=Any)
    model = derive_single_flavor_gnw_tensor()
    rules = _substitution_rules(model, parameters, substitutions)
    return _coefficient_array_result(model.tensor, _pair_groups(model.legs);
        materialize, substitutions=rules, simplify, element_type)
end

function schwinger_fermion_coefficient_tensor(; parameters=(;), substitutions=nothing,
        materialize::Bool=true, simplify::Bool=false, element_type=Any)
    model = derive_schwinger_fermion_tensor()
    rules = _substitution_rules(model, parameters, substitutions)
    return _coefficient_array_result(model.tensor, _pair_groups(model.legs);
        materialize, substitutions=rules, simplify, element_type)
end

function two_color_qcd_fermion_coefficient_tensor(; parameters=(;), substitutions=nothing,
        diquark_source::Bool=false, materialize::Bool=true, simplify::Bool=false, element_type=Any)
    model = derive_two_color_qcd_fermion_tensor(; diquark_source)
    rules = _substitution_rules(model.vars, parameters, substitutions)
    return _coefficient_array_result(model.tensor, _four_groups_16(model.legs);
        materialize, substitutions=rules, simplify, element_type)
end

function two_flavor_gnw_leg_groups(model)
    return _four_groups_16(model.legs)
end

function two_flavor_gnw_physical_substitutions(model; M1, M2, gs, gpi,
        r1=1, r2=1, mu1=0, mu2=0, H=0, q=1 / sqrt(2))
    # Ref/2d_gn.pdf section 3 sets Wilson r^(f)=1 in the hopping projectors.
    # For r != 1 this convenience mapping only shifts the onsite Wilson mass
    # M_f+2 -> M_f+2r_f; the hopping decomposition itself remains the r=1 one.
    p = model.params
    return Dict(
        p.M1 => M1 + 2r1 - 2,
        p.M2 => M2 + 2r2 - 2,
        p.gs2 => gs^2,
        p.gp2 => gpi^2,
        p.H => H,
        p.q => q,
        p.z1 => exp(mu1 / 2),
        p.z2 => exp(mu2 / 2),
    )
end

function two_flavor_gnw_coefficient_tensors(; parameters=(;), physical_parameters=nothing,
        substitutions=nothing, materialize::Bool=true, simplify::Bool=false, element_type=Any)
    names = (:pure, :chiral_condensate, :pseudoscalar_singlet, :pseudoscalar_triplet)

    if !materialize
        @variables M1 M2 gs2 gp2 H q z1 z2
        params = (; M1, M2, gs2, gp2, H, q, z1, z2)
        vars = _tfgnw_generators()
        lightweight_model = (; params, legs=vars.legs)
        groups = two_flavor_gnw_leg_groups(lightweight_model)
        rules = substitutions
        if rules === nothing
            rules = physical_parameters === nothing ?
                _substitution_rules(lightweight_model, parameters, nothing) :
                two_flavor_gnw_physical_substitutions(lightweight_model; physical_parameters...)
        end
        values = map(names) do name
            _coefficient_array_result(name, groups;
                materialize=false, substitutions=rules, simplify, element_type)
        end
        return NamedTuple{names}(values)
    end

    model = derive_two_flavor_gnw_tensors()
    groups = two_flavor_gnw_leg_groups(model)
    rules = substitutions
    if rules === nothing
        rules = physical_parameters === nothing ?
            _substitution_rules(model, parameters, nothing) :
            two_flavor_gnw_physical_substitutions(model; physical_parameters...)
    end
    values = map(names) do name
        _coefficient_array_result(getproperty(model, name), groups;
            materialize, substitutions=rules, simplify, element_type)
    end
    return NamedTuple{names}(values)
end

function three_flavor_gnw_leg_groups(model)
    return [collect(model.legs[1:4]), collect(model.legs[5:8]), collect(model.legs[9:12]),
        collect(model.legs[13:16]), collect(model.legs[17:20]), collect(model.legs[21:24])]
end

function three_flavor_gnw_coefficient_tensors(; parameters=(;), substitutions=nothing,
        observables=(:pure, :pseudoscalar_singlet_squared, :pseudoscalar_lambda8_squared),
        materialize::Bool=false, simplify::Bool=false, element_type=Any)
    model = three_flavor_gnw_model()
    groups = three_flavor_gnw_leg_groups(model)
    rules = _substitution_rules(model, parameters, substitutions)
    values = map(observables) do observable
        oracle = bits -> three_flavor_gnw_coefficient(model, bits;
            observable=observable,
            simplify=false,
        )
        _coefficient_array_result(oracle, groups;
            materialize, substitutions=rules, simplify, element_type)
    end
    return NamedTuple{Tuple(observables)}(values)
end

function njl_coefficient_tensor(; parameters=(;), substitutions=nothing,
        n1=0, n2=0, n3=0, materialize::Bool=false, simplify::Bool=false, element_type=Any)
    model = njl_variables()
    p = model.parameters
    rules = _substitution_rules(model, parameters, substitutions)
    oracle = bits_tuple -> begin
        bits = njl_bits(bits_tuple)
        njl_I_coefficient(bits, p.m, p.r, p.u, n1, n2, n3) * njl_S_sign(bits)
    end
    return _coefficient_array_result(oracle, (3, 3, 3, 3, 3, 3, 3, 3);
        materialize, substitutions=rules, simplify, element_type)
end

function hubbard_coefficient_tensor(; parameters=(;), substitutions=nothing,
        materialize::Bool=true, simplify::Bool=false, element_type=Any)
    model = derive_hubbard_tensor()
    groups = [collect(model.legs[1:4]), collect(model.legs[5:6]),
        collect(model.legs[7:8]), collect(model.legs[9:12])]
    rules = _substitution_rules(model, parameters, substitutions)
    return _coefficient_array_result(model.tensor, groups;
        materialize, substitutions=rules, simplify, element_type)
end
