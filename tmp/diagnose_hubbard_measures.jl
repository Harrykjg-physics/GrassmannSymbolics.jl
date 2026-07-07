using Pkg

slash = Char(47)
dot = Char(46)
root = @__DIR__() * string(slash, dot, dot)
Pkg.activate(root)

include(joinpath(root, string(:examples), string(1, :D_Hubbard, dot, :jl)))

base = hubbard_variables()
p = base.parameters
o = base.original
l = base.named_legs

function tensor_with_measure(order)
    nup = o.ψbaru * o.ψu
    ndn = o.ψbard * o.ψd
    onsite_action = p.Ue * nup * ndn - p.h * (nup + ndn)
    bond_exponents = [
        p.q * o.ψbaru * l.ησu,
        p.q * o.ψbard * l.ησd,
        p.q * l.ζσu * o.ψu,
        p.q * l.ζσd * o.ψd,
        -o.ψbaru * l.ητu,
        -o.ψbard * l.ητd,
        l.ηbarτu * o.ψu,
        l.ηbarτd * o.ψd,
        -p.q * o.ψbaru * l.ζbarσu,
        -p.q * o.ψbard * l.ζbarσd,
        p.q * l.ηbarσu * o.ψu,
        p.q * l.ηbarσd * o.ψd,
    ]
    integrand = exp(-onsite_action)
    for term in bond_exponents
        integrand = integrand * exp(term)
    end
    return integrate(integrand, order)
end

orders = (
    (o.ψbaru, o.ψu, o.ψbard, o.ψd),
    (o.ψu, o.ψbaru, o.ψd, o.ψbard),
    (o.ψbard, o.ψd, o.ψbaru, o.ψu),
    (o.ψd, o.ψbard, o.ψu, o.ψbaru),
    (o.ψbaru, o.ψbard, o.ψu, o.ψd),
    (o.ψu, o.ψd, o.ψbaru, o.ψbard),
)

for order in orders
    tensor = tensor_with_measure(collect(order))
    mismatches = 0
    first = nothing
    for mask in 0:4095
        leg_bits = hubbard_leg_bits(mask)
        bits = hubbard_bits_from_leg_bits(leg_bits)
        occupied = GrassmannGenerator[
            leg for (leg, occupation) in zip(base.legs, leg_bits) if occupation == 1
        ]
        actual = get_coeff(tensor, occupied)
        expected = published_hubbard_coefficient(bits, p.q, p.h, p.Ue, base.named_legs)
        if !symbolically_equal(GrassmannExpr(actual), GrassmannExpr(expected))
            mismatches += 1
            first === nothing && (first = (mask, actual, expected))
        end
    end
    println(string(:order, :_, order))
    println(string(:mismatches, :_, mismatches))
    println(string(:first, :_, first))
end
