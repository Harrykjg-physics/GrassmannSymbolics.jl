using Pkg

slash = Char(47)
dot = Char(46)
root = @__DIR__() * string(slash, dot, dot)
Pkg.activate(root)

include(joinpath(root, string(:examples), string(:NJL, dot, :jl)))

model = njl_variables()
p = model.parameters
candidates = njl_nonzero_candidate_bits()

for n1 in 0:1, n2 in 0:1, n3 in 0:1
    for bits in candidates
        direct = njl_direct_local_integral(bits, p.m, p.r, p.u, n1, n2, n3) *
            njl_raw_G_monomial(bits, model.legs)
        published = njl_projected_entry_without_measure_sign(bits, model; n1=n1, n2=n2, n3=n3)
        if !symbolically_equal(direct, published)
            println(string(:mismatch))
            println(bits)
            println(string(:parity, :_, n1, :_, n2, :_, n3))
            println(string(:bar_degree, :_, njl_bar_degree(bits)))
            println(string(:chi_degree, :_, njl_chi_degree(bits)))
            println(string(:S, :_, njl_S_sign(bits)))
            println(string(:I, :_, njl_I_coefficient(bits, p.m, p.r, p.u, n1, n2, n3)))
            println(string(:direct, :_, direct))
            println(string(:published, :_, published))
            println(string(:difference, :_, simplify_coefficients(direct - published)))
            exit()
        end
    end
end

println(string(:no_mismatch))
