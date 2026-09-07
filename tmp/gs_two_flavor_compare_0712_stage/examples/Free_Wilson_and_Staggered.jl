"""
Symbolic derivation of the fundamental Grassmann tensors in appendix B of
Kadoh et al., JHEP 10 (2021) 188.

The implementation follows equations (3.7), (3.8), and (3.10). It keeps the
Berezin measure and tensor-leg ordering explicit and audits every coefficient
against the published Wilson tensor (B.4) and staggered tensor (B.6). The
published expressions are retained separately when the audit detects a mismatch.
"""

using GrassmannSymbolics
using Symbolics
using Test

function linear_combination(coefficients, fields)
    result = GrassmannExpr(0)
    for (coefficient, field) in zip(coefficients, fields)
        result = result + coefficient * field
    end
    return result
end

"""Return the local exponent factors in equation (3.10) for one direction."""
function directional_bond_exponents(
    ψ,
    ψbar,
    U,
    Vt,
    η,
    ζ,
    ηbar,
    ζbar;
    x_channel=2,
    y_channel=1,
    σx=1,
    σy=1,
)
    x_left = linear_combination(U[:, x_channel], ψbar)
    x_right = linear_combination(Vt[x_channel, :], ψ)
    y_left = linear_combination(U[:, y_channel], ψbar)
    y_right = linear_combination(Vt[y_channel, :], ψ)

    x_exponent = -x_left * η + σx * ηbar * x_right
    y_exponent = y_left * ζbar + σy * ζ * y_right
    return (x_exponent, y_exponent)
end

function derive_wilson_tensor(; x_channel=2, y_channel=1)
    @variables D A11 A12 A21 A22 B11 B22

    ψ = grassmann(:ψ, 1:2)
    ψbar = dual.(ψ)
    ηx, ζx, ηt, ζt = (
        grassmann(:ηx), grassmann(:ζx), grassmann(:ηt), grassmann(:ζt)
    )
    ηbarx, ζbarx, ηbart, ζbart = dual.((ηx, ζx, ηt, ζt))

    # A = U_x, B = U_t and V_μ† = -U_μ in equations (B.2)-(B.3).
    # The paper's printed (B.4) uses the equivalent temporal SVD phase
    # U_t -> -U_t, V_t† -> -V_t†, i.e. ηt,ζt -> -ηt,-ζt.
    A = [A11 A12; A21 A22]
    B = [B11 0; 0 B22]
    onsite_action = D * (ψbar[1] * ψ[1] + ψbar[2] * ψ[2])
    bond_exponents = [
        directional_bond_exponents(
            ψ, ψbar, A, -A, ηx, ζx, ηbarx, ζbarx;
            x_channel, y_channel,
        )...,
        directional_bond_exponents(
            ψ, ψbar, -B, B, ηt, ζt, ηbart, ζbart;
            x_channel, y_channel,
        )...,
    ]

    measure_order = [ψ[1], ψbar[1], ψ[2], ψbar[2]]
    tensor = local_grassmann_tensor(
        onsite_action,
        bond_exponents,
        measure_order,
    )

    # Equation (B.4), transcribed without algebraic rearrangement.
    first_bracket = (
        ζt * ηbart - ζx * ηbarx +
        A11 * B22 * ζt * ηbarx +
        B11 * A21 * ηbart * ηbarx +
        A12 * B22 * ζx * ζt +
        B11 * A22 * ζx * ηbart
    )
    second_bracket = (
        ηt * ζbart - ηx * ζbarx +
        B22 * A11 * ηx * ζbart -
        A12 * B11 * ηx * ηt -
        B22 * A21 * ζbart * ζbarx +
        A22 * B11 * ηt * ζbarx
    )
    reference = (
        D^2 -
        D * B11^2 * ηt * ηbart -
        D * B22^2 * ζt * ζbart -
        D * (A21 * A12 + A11^2) * ηx * ηbarx -
        D * (A22^2 + A12 * A21) * ζx * ζbarx +
        D * (A22 * A12 + A12 * A11) * ηx * ζx -
        D * (A21 * A22 + A11 * A21) * ζbarx * ηbarx +
        D * B22 * A12 * ηx * ζt -
        D * A21 * B22 * ζbart * ηbarx -
        D * A22 * B22 * ζx * ζbart -
        D * B22 * A22 * ζt * ζbarx -
        D * A12 * B11 * ζx * ηt +
        D * B11 * A21 * ηbart * ζbarx -
        D * B11 * A11 * ηx * ηbart -
        D * A11 * B11 * ηt * ηbarx -
        first_bracket * second_bracket
    )

    legs = (
        Ψx=(ηx, ζx),
        Ψt=(ηt, ζt),
        Ψbar_t=(ζbart, ηbart),
        Ψbar_x=(ζbarx, ηbarx),
    )
    parameters = (; D, A11, A12, A21, A22, B11, B22)
    return (; tensor, published_tensor=simplify_coefficients(reference), legs, parameters)
end

function derive_staggered_tensor()
    @variables M p

    χ = grassmann(:χ)
    χbar = dual(χ)
    ηx, ζx, ηt, ζt = (
        grassmann(:ηx), grassmann(:ζx), grassmann(:ηt), grassmann(:ζt)
    )
    ηbarx, ζbarx, ηbart, ζbart = dual.((ηx, ζx, ηt, ζt))

    # X_μ=p_μ/2 and Y_μ=-p_μ/2. We place the staggered phase on
    # the barred-field factor of X and on the unbarred-field factor of Y.
    x_exponents = (
        -χbar * ηx + (1 // 2) * ηbarx * χ,
        χbar * ζbarx - (1 // 2) * ζx * χ,
    )
    t_exponents = (
        -p * χbar * ηt + (1 // 2) * ηbart * χ,
        χbar * ζbart - (p * (1 // 2)) * ζt * χ,
    )
    tensor = local_grassmann_tensor(
        M * χbar * χ,
        [x_exponents..., t_exponents...],
        [χ, χbar],
    )

    # Equation (B.6), with p=p_t(n)=(-1)^n_x.
    reference = (
        -M -
        (1 // 2) * ηx * ηbarx +
        (1 // 2) * ζx * ζbarx -
        (p * (1 // 2)) * ηt * ηbart +
        (p * (1 // 2)) * ζt * ζbart -
        (1 // 2) * ηx * ζx -
        (1 // 2) * ηt * ζt -
        (1 // 2) * ζbarx * ηbarx -
        (1 // 2) * ζbart * ηbart -
        (p * (1 // 2)) * ηx * ζt +
        (p * (1 // 2)) * ζx * ηt +
        (1 // 2) * ζbart * ηbarx -
        (1 // 2) * ηbart * ζbarx -
        (1 // 2) * ηx * ηbart +
        (1 // 2) * ζx * ζbart -
        (p * (1 // 2)) * ηt * ηbarx +
        (p * (1 // 2)) * ζt * ζbarx
    )

    legs = (
        Ψx=(ηx, ζx),
        Ψt=(ηt, ζt),
        Ψbar_t=(ζbart, ηbart),
        Ψbar_x=(ζbarx, ηbarx),
    )
    return (; tensor, published_tensor=simplify_coefficients(reference), legs)
end

"""Close all four legs of a one-site tensor with the paper's Gaussian measure."""
function contract_periodic_one_site(tensor, legs)
    ηx, ζx = legs.Ψx
    ηt, ζt = legs.Ψt
    ζbart, ηbart = legs.Ψbar_t
    ζbarx, ηbarx = legs.Ψbar_x

    pairs = (
        (ηbarx, ηx),
        (ζbarx, ζx),
        (ηbart, ηt),
        (ζbart, ζt),
    )
    return scalar_part(contract_grassmann((tensor,), pairs))
end

function specialize_staggered_parity(expression, p, parity)
    return map_coefficients(expression) do coefficient
        try
            Symbolics.substitute(coefficient, Dict(p => parity))
        catch
            coefficient
        end
    end
end

function maximum_abs_numeric_coefficient(expression)
    maximum_value = BigFloat(0)
    for coefficient in values(expression.terms)
        value = try
            Symbolics.value(coefficient)
        catch
            coefficient
        end
        maximum_value = max(maximum_value, abs(BigFloat(value)))
    end
    return maximum_value
end

function wilson_literature_error(wilson)
    return setprecision(256) do
        q = inv(sqrt(BigFloat(2)))
        p = wilson.parameters
        maximum_error = BigFloat(0)
        for onsite_value in (BigFloat(-5) / 4, BigFloat(0), BigFloat(11) / 8)
            rules = Dict(
                p.D => onsite_value,
                p.A11 => -q,
                p.A12 => -q,
                p.A21 => -q,
                p.A22 => q,
                p.B11 => BigFloat(-1),
                p.B22 => BigFloat(-1),
            )
            difference = simplify_coefficients(substitute_coefficients(
                wilson.tensor - wilson.published_tensor,
                rules,
            ))
            maximum_error = max(
                maximum_error,
                maximum_abs_numeric_coefficient(difference),
            )
        end
        maximum_error
    end
end
function validate_free_tensors()
    @testset "Free Wilson and staggered fundamental tensors" begin
        wilson = derive_wilson_tensor()
        staggered = derive_staggered_tensor()

        # Equation (B.4) is written after substituting the fixed matrices in
        # (B.3), so compare at 256-bit precision for three independent D values.
        @test wilson_literature_error(wilson) < big"1e-60"

        @variables p
        for parity in (-1, 1)
            tensor_at_parity = specialize_staggered_parity(
                staggered.tensor, p, parity,
            )
            reference_at_parity = specialize_staggered_parity(
                staggered.published_tensor, p, parity,
            )
            parity_ok = symbolically_equal(tensor_at_parity, reference_at_parity)
            parity_ok || println(
                "Staggered difference at p=", parity, ": ",
                simplify_coefficients(tensor_at_parity - reference_at_parity),
            )
            @test_broken parity_ok

            # A one-component local integral has quadratic part A*B and is
            # therefore a decomposable two-form. Printed (B.6) fails this
            # necessary condition as well as the coefficient comparison.
            reference_quadratic = reference_at_parity - scalar_part(reference_at_parity)
            @test_broken iszero(reference_quadratic^2)
        end
        @test is_even(wilson.tensor)
        @test is_even(staggered.tensor)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    validate_free_tensors()
    println("Wilson tensor:   ", derive_wilson_tensor().tensor)
    println("Staggered tensor: ", derive_staggered_tensor().tensor)
end
