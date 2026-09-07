
########################### Grassmann tensors for (1+1)D lattice models ###########################

"""
Grassmann (coefficient) tensor for (1+1)D spinless free fermion

    H = -t ∑(ci†cj + h.c.) - μ ∑ni

Arguments: 
    t : hopping parameter
    μ : chemical potential
    ϵ : discretization in the temporal direction
Return :
    T_coef : A rank-4 bulk tensor : T_coef[(iσ, jσ), iτ, (iσp, jσp), iτp]
    T_coef_n_imp : A rank-4 n impurity tensor  
    T_coef_cdag_imp : A rank-4 c+ impurity tensor 
    T_coef_c_imp : A rank-4 c impurity tensor 

    Notice that the order of Grassmann numbers is not consistent with the order of coefficient indices
    Therefore add_perm_sign(T_coef, (1, 2, 4, 3)) is needed for further auto_sign functions
"""

function free_spinless_fermion_1d(t, μ, ϵ)

    Q = eltype(t)
    T_coef = zeros(Q, (4, 2, 4, 2))
    T_coef_n_imp = zeros(Q, (4, 2, 4, 2))
    T_coef_cdag_imp = zeros(Q, (4, 2, 4, 2))
    T_coef_c_imp = zeros(Q, (4, 2, 4, 2))

    for iσ in 0:1, jσ in 0:1, iτ in 0:1, iσp in 0:1, jσp in 0:1, iτp in 0:1

        Iσ = f(iσ, jσ); Iσp = f(iσp, jσp)
        Iτ = f(iτ); Iτp = f(iτp)

        exponent1 = iσp + jσ + jσp + iτ + iτp
        exponent2 = iσ + iσp + jσ + jσp
        exponent3 = jσp * (iσ + iτp + jσ) + iτp * jσ + iτ * (iσ + jσ)
        sign_num = (-1)^exponent3 
        # sign_num = auto_sign((iτ, jσp, iσ, iτp, jσ, iσp), (3, 5, 1, 4, 2, 6))

        T_coef[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num *
                                                    (isequal(iσp + jσ + iτp, 1) * isequal(iσ + jσp + iτ, 1) -
                                                    ((μ * ϵ + 1) * isequal(iσp + jσ + iτp, 0) * isequal(iσ + jσp + iτ, 0)))

        T_coef_n_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num *
                                            (isequal(iσp + jσ + iτp, 0) * isequal(iσ + jσp + iτ, 0))

        T_coef_cdag_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num *
                                             isequal(iσp + jσ + iτp, 1) * isequal(iσ + jσp + iτ, 0)

        T_coef_c_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num *
                                             isequal(iσp + jσ + iτp, 0) * isequal(iσ + jσp + iτ, 1)
    end

    return T_coef, T_coef_n_imp, T_coef_cdag_imp, T_coef_c_imp
end

# ---------------------------------------------------------------------------

"""
Grassmann (coefficient) tensor for (1+1)D Hubbard model

    H = -t ∑(ci†cj + h.c.) - μ ∑ni + U/2 ∑ ni↑ni↓

Arguments: 
    t : hopping parameter
    μ : chemical potential
    U : Hubbard interaction
    ϵ : discretization in the temporal direction

Return :
    T_coef : A rank-4 bulk tensor : T_coef[(iσu, iσd, jσu, jσd), (iτu, iτd), (iσup, iσdp, jσup, jσdp), (iτup, iτdp)]
    T_coef_n_imp : A rank-4 n impurity tensor 
    T_coef_cdagup_imp : A rank-4 cdagup impurity tensor 
    T_coef_cup_imp : A rank-4 cup impurity tensor 
    T_coef_cdagdn_imp : A rank-4 cdagdn impurity tensor 
    T_coef_cdn_imp : A rank-4 cdn impurity tensor 

    Notice that the order of Grassmann numbers is not consistent with the order of coefficient indices
    Therefore add_perm_sign(T_coef, (1, 2, 4, 3)) is needed for further auto_sign functions
"""

function Hubbard_model_1d(t, μ, ϵ, U)

    Q = eltype(t)
    T_coef = zeros(Q, (16, 4, 16, 4))
    T_coef_n_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cdagup_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cup_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cdagdn_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cdn_imp = zeros(Q, (16, 4, 16, 4))

    for iσu in 0:1, iσd in 0:1, jσu in 0:1, jσd in 0:1, iτu in 0:1, iτd in 0:1
        for iσup in 0:1, iσdp in 0:1, jσup in 0:1, jσdp in 0:1, iτup in 0:1, iτdp in 0:1

            Iσ = f(iσu, iσd, jσu, jσd); Iσp = f(iσup, iσdp, jσup, jσdp)
            Iτ = f(iτu, iτd); Iτp = f(iτup, iτdp)

            exponent1 = iσu + iσd
            exponent2 = iσu + jσu + iσup + jσup + iσd + jσd + iσdp + jσdp
            exponent3 = iσu * iτu + iσd * (iτu + jσup + iτup + iσup + jσu + iτd) + jσu * (iτu + jσup + iτup + iσup) +
                        jσd * (iτu + jσup + iτup + iσup + iτd + jσdp + iτdp + iσdp) + iτd * (jσup + iτup + iσup) +
                        iτdp * (jσup + iτup + iσup + jσdp) + iτup * jσup +
                        jσdp * (jσup + iσup) + iσdp * iσup
            sign_num = (-1)^exponent3
            # sign_num = auto_sign((iτu, iσu, jσup, iτup, iσup, jσu, iτd, iσd, jσdp, iτdp, iσdp, jσd), (2, 8, 6, 12, 1, 7, 10, 4, 9, 3, 11, 5))

            T_coef[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
                                   (isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1) -
                                   (μ * ϵ + 1) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1)) -
                                   (μ * ϵ + 1) * (isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0)) -
                                   (U * ϵ - (μ * ϵ + 1)^2) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0))
                                    )

            T_coef_n_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
                                       (
                                       -(isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0)) -
                                       (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1)) +
                                       2 * (μ * ϵ + 1) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0))
                                       ) 

            T_coef_cdagup_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
                                                (isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 1) -
                                                (μ * ϵ + 1) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 1))
                                                )

            T_coef_cup_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
                                                (isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 0) -
                                                (μ * ϵ + 1) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 0))
                                                )

            T_coef_cdagdn_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
                                                (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1) -
                                                (μ * ϵ + 1) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0))
                                                )

            T_coef_cdn_imp[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
                                                (isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1) -
                                                (μ * ϵ + 1) * (isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0))
                                                )

        end
    end

    return T_coef, T_coef_n_imp, T_coef_cdagup_imp, T_coef_cup_imp, T_coef_cdagdn_imp, T_coef_cdn_imp
end

# ---------------------------------------------------------------------------

"""
Grassmann (coefficient) tensor for (1+1)D Hubbard model with chemical potential and B field

    H = -t ∑(cis†cjs + h.c.) + U ∑ (ni↑-1/2)(ni↓-1/2) - μ(ni↑ + ni↓) - B(ni↑ - ni↓)

Arguments: 
    t : hopping parameter
    μ : chemical potential
    B : applied magnetic field
    U : Hubbard interaction
    ϵ : discretization in the temporal direction
"""

function Hubbard_model_1d(t, μ, B, ϵ, U)

    Q = eltype(t)
    T_coef = zeros(Q, (16, 4, 16, 4))
    T_coef_nu_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_nd_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_n_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cdagup_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cdagdn_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cup_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_cdn_imp = zeros(Q, (16, 4, 16, 4))
    T_coef_nupndn_imp = zeros(Q, (16, 4, 16, 4))

    for iσu in 0:1, iσd in 0:1, jσu in 0:1, jσd in 0:1, iτu in 0:1, iτd in 0:1
        for iσup in 0:1, iσdp in 0:1, jσup in 0:1, jσdp in 0:1, iτup in 0:1, iτdp in 0:1

            Iσ = f(iσu, iσd, jσu, jσd); Iσp = f(iσup, iσdp, jσup, jσdp)
            Iτ = f(iτu, iτd); Iτp = f(iτup, iτdp)

            exponent1 = iσu + iσd
            exponent2 = iσu + jσu + iσup + jσup + iσd + jσd + iσdp + jσdp
            exponent3 = iσu * jσup + iσd * (jσdp + iτup + iσup + jσu + iτu + jσup) + jσu * (iτu + jσup) + jσd * (iτd + jσdp + iτup + iσup + iτu + jσup)
            + iτu * jσup + iτd * (jσdp + iτup + iσup + jσup) + iτdp * (iσdp + jσdp + iτup + iσup + jσup) + iτup * (iσup + jσup) +
            jσdp * (iσup + jσup) + iσdp * iσup

            # sign_num = (-1)^exponent3
            sign_num = auto_sign((jσup, iσu, iτu, jσu, iσup, iτup, jσdp, iσd, iτd, jσd, iσdp, iτdp), (2, 8, 4, 10, 3, 9, 12, 6, 7, 1, 11, 5))
            T_coef[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            (
                isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1) -
                (ϵ * (μ + B + U/2) + 1) * (isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0)) -
                (ϵ * (μ - B + U/2) + 1) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1)) -
                (ϵ * U - (ϵ * (μ + B + U/2) + 1) * (ϵ * (μ - B + U/2) + 1)) * (isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0))
            )

            T_coef_nu_imp[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            (
                - isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0)  +
                (ϵ * (μ - B + U/2) + 1) * isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0) 
            )

            T_coef_nd_imp[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            (
                - isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1)  +
                (ϵ * (μ + B + U/2) + 1) * isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0) 
            )

            T_coef_n_imp = T_coef_nu_imp + T_coef_nd_imp

            T_coef_cdagup_imp[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            (
                isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 1)  +
                (ϵ * (μ - B + U/2) + 1) * isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 1) 
            )

            T_coef_cdagdn_imp[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            (
                isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1)  +
                (ϵ * (μ + B + U/2) + 1) * isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0) 
            )

            T_coef_cup_imp[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            (
                - isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 1) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 0)  +
                (ϵ * (μ - B + U/2) + 1) * isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 0) 
            )

            T_coef_cdn_imp[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            (
                - isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 1) * isequal(iτup + iσup + jσu, 1)  +
                (ϵ * (μ + B + U/2) + 1) * isequal(iτd + iσd + jσdp, 1) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0) 
            )

            T_coef_nupndn_imp[Iσ, Iτ, Iσp, Iτp] = exp(-U * ϵ / 4) * (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * 
            isequal(iτd + iσd + jσdp, 0) * isequal(iτdp + iσdp + jσd, 0) * isequal(iτu + iσu + jσup, 0) * isequal(iτup + iσup + jσu, 0)
            
        end
    end

    return T_coef, T_coef_nu_imp, T_coef_nd_imp, T_coef_n_imp, T_coef_cdagup_imp, T_coef_cdagdn_imp, T_coef_cup_imp, T_coef_cdn_imp, T_coef_nupndn_imp
end

# ---------------------------------------------------------------------------

"""
Grassmann (coefficient) tensor for (1+1)D extended Hubbard model

    H = -t ∑(ci†cj + h.c.) + U ∑ (ni↑-1/2)(ni↓-1/2) + V ∑ (ni-1)(nj-1) - μ ∑((ni↑ + ni↓) + B(ni↑ - ni↓)) 

Arguments: 
    t : hopping paramter
    μ : chemical potential
    U : on-site Hubbard interaction
    V : nearest-neighbour Hubbard interaction
    B : external magnetic field
    ϵ : discretization in the temporal direction
Return :
    T_coef : A rank-4 bulk tensor : T_coef[(iσu, iσd, jσu, jσd), (iτu, iτd), (iσup, iσdp, jσup, jσdp), (iτup, iτdp)]

    Notice that the order of Grassmann numbers is not consistent with the order of coefficient indices
    Therefore add_perm_sign(T_coef, (1, 2, 4, 3)) is needed for further auto_sign functions
"""

function extended_Hubbard_model_1d(t, μ, ϵ, U, V, B)

    N1 = -1 - (U*ϵ)/2 - V*ϵ - μ*ϵ - B*ϵ
    N2 = -1 - (U*ϵ)/2 - V*ϵ - μ*ϵ + B*ϵ
    N3 = (U*ϵ)/4 + V*ϵ + μ*ϵ

    # T_coef[(iσu, iσd, jσu, jσd, kσu), (iτu, iτd), (iσup, iσdp, jσup, jσdp, kσup), (iτup, iτdp)]
    T_coef = zeros(eltype(t), (32, 4, 32, 4))

    for iτdp in 0:1, iτup in 0:1, kσup in 0:1, jσdp in 0:1, jσup in 0:1, iσdp in 0:1, iσup in 0:1
        for iτd in 0:1, iτu in 0:1, kσu in 0:1, jσd in 0:1, jσu in 0:1, iσd in 0:1, iσu in 0:1

            Iσ = f(iσu, iσd, jσu, jσd, kσu)
            Iτ = f(iτu, iτd)
            Iσp = f(iσup, iσdp, jσup, jσdp, kσup)
            Iτp = f(iτup, iτdp)

            exponent1 = iτu + iτd + iσup + iσdp + iτup + iτdp + jσu + jσd + jσup + jσdp
            exponent2 = iσu + jσu + iσup + jσup + iσd + jσd + iσdp + jσdp
            exponent3 = kσu + kσup

            sign_perm = auto_sign((iσu, iσd, jσu, jσd, kσu, iτu, iτd, iσup, iσdp, jσup, jσdp, kσup, iτup, iτdp), (2, 8, 4, 10, 14, 3, 9, 13, 7, 1, 11, 5, 12, 6))

            T_coef[Iσ, Iτ, Iσp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * (sqrt(-ϵ * V))^exponent3 * exp(-N3) * (
            isequal(iτdp + iσdp + jσd, 1) *  isequal(iτd + iσd + jσdp, 1) * isequal(iτup + iσup + jσu + kσu + kσup, 1) * isequal(iτu + iσu + jσup + kσu + kσup, 1)
            - N1 * (isequal(iτdp + iσdp + jσd, 1) *  isequal(iτd + iσd + jσdp, 1) * isequal(iτup + iσup + jσu + kσu + kσup, 0) * isequal(iτu + iσu + jσup + kσu + kσup, 0))
            - N2 * (isequal(iτdp + iσdp + jσd, 0) *  isequal(iτd + iσd + jσdp, 0) * isequal(iτup + iσup + jσu + kσu + kσup, 1) * isequal(iτu + iσu + jσup + kσu + kσup, 1))
            + (N1*N2 - U*ϵ) * (isequal(iτdp + iσdp + jσd, 0) *  isequal(iτd + iσd + jσdp, 0) * isequal(iτup + iσup + jσu + kσu + kσup, 0) * isequal(iτu + iσu + jσup + kσu + kσup, 0))
            ) * sign_perm
                                                                                                           
        end
    end

    return T_coef
end

# ---------------------------------------------------------------------------

"""
Grassmann (coefficient) tensor for the (1+1)D SSH model with type-II interaction

    H = -J1 ∑(ciA†cjB + h.c.) - J2 ∑(ciB†cjA + h.c.) + U ∑ ni nj

Arguments: 
    J1 : hopping paramter within a unit cell
    J2 : hopping paramter between unit cells
    V : nearest-neighbour interaction
    ϵ : discretization in the temporal direction
Return :
    T_coef_A : A rank-4 bulk tensor : T_coef_A[(iσA, jσA, kσA), iτA, (iσBp, jσBp, kσBp), iτAp]
    T_coef_B : A rank-4 bulk tensor : T_coef_B[(iσB, jσB, kσB), iτB, (iσAp, jσAp, kσAp), iτBp]
"""

function SSH_model_typeII(J1, J2, V, ϵ)

    # T_coef_A[(iσA, jσA, kσA), iτA, (iσBp, jσBp, kσBp), iτAp]
    T_coef_A = zeros(eltype(J1), (8, 2, 8, 2))
    # T_coef_B[(iσB, jσB, kσB), iτB, (iσAp, jσAp, kσAp), iτBp]
    T_coef_B = zeros(eltype(J1), (8, 2, 8, 2))

    for iσA in 0:1, jσA in 0:1, kσA in 0:1, iτA in 0:1
        for iσBp in 0:1, jσBp in 0:1, kσBp in 0:1, iτAp in 0:1

            IσA = f(iσA, jσA, kσA)
            IτA = f(iτA)
            IσBp = f(iσBp, jσBp, kσBp)
            IτAp = f(iτAp)

            exponent1 = iσA + kσBp * kσA
            exponent2 = iσA + jσA
            exponent3 = iσBp + jσBp
            exponent4 = kσA + kσBp
            S = iσA * iτA + jσA * (kσA + kσBp + jσBp + iτA) + kσA * (kσBp + jσBp + iτA) + kσBp * jσBp
            sign_perm = (-1)^S
            # sign_perm = auto_sign((iτA, iσA, jσBp, kσBp, kσA, jσA, iσBp, iτAp), (2, 6, 5, 1, 4, 3, 7, 8))
            T_coef_A[IσA, IτA, IσBp, IτAp] = (-1)^exponent1 * (sqrt(ϵ * J1))^exponent2 * (sqrt(ϵ * J2))^exponent3 * (sqrt(ϵ * V))^exponent4 * (
            isequal(iτAp + iσBp + jσA + kσA + kσBp, 1) *  isequal(kσA + kσBp + jσBp + iσA + iτA, 1) - 
            isequal(iτAp + iσBp + jσA + kσA + kσBp, 0) *  isequal(kσA + kσBp + jσBp + iσA + iτA, 0)
            ) * sign_perm
                                                                                                           
        end
    end

    for iσB in 0:1, jσB in 0:1, kσB in 0:1, iτB in 0:1
        for iσAp in 0:1, jσAp in 0:1, kσAp in 0:1, iτBp in 0:1

            IσB = f(iσB, jσB, kσB)
            IτB = f(iτB)
            IσAp = f(iσAp, jσAp, kσAp)
            IτBp = f(iτBp)

            exponent1 = iσB
            exponent2 = iσAp + jσAp
            exponent3 = iσB + jσB
            exponent4 = kσAp + kσB
            S = iσB * (jσAp + iτB) + jσB * (kσB + kσAp + jσAp + iτB) + kσB * (kσAp + jσAp + iτB) + kσAp * jσAp
            sign_perm = (-1)^S
            # sign_perm = auto_sign((iτB, jσAp, iσB, kσAp, kσB, jσB, iσAp, iτBp), (3, 6, 5, 1, 4, 2, 7, 8))
            T_coef_B[IσB, IτB, IσAp, IτBp] = (-1)^exponent1 * (sqrt(ϵ * J1))^exponent2 * (sqrt(ϵ * J2))^exponent3 * (sqrt(ϵ * V))^exponent4 * (
            isequal(iτBp + iσAp + jσB + kσAp + kσB, 1) *  isequal(kσB + kσAp + jσAp + iσB + iτB, 1) -
            isequal(iτBp + iσAp + jσB + kσAp + kσB, 0) *  isequal(kσB + kσAp + jσAp + iσB + iτB, 0)
            ) * sign_perm
                                                                                                           
        end
    end

    return T_coef_A, T_coef_B
end

########################### Grassmann tensors for (2+1)D lattice models ###########################

"""
Grassmann (coefficient) tensor for (2+1)D Hubbard model on the square lattice

    H = -t ∑(ci†cj + h.c.) - μ ∑ni + U/2 ∑ ni↑ni↓

Arguments: 
    t : hopping paramter
    μ : chemical potential
    U : Hubbard interaction
    ϵ : discretization in the temporal direction
Return :
    T_coef : A rank-6 bulk tensor : T_coef[(ixu, ixd, jxu, jxd), (iyu, iyd, jyu, jyd), (iτu, iτd), (ixup, ixdp, jxup, jxdp), (iyup, iydp, jyup, jydp), (iτup, iτdp)]

    Notice that the order of Grassmann numbers is not consistent with the order of coefficient indices
    Therefore add_perm_sign(T_coef, (1, 2, 3, 6, 5, 4)) is needed for further auto_sign functions
"""

function Hubbard_model_square(t, μ, ϵ, U)

    T_coef = zeros(eltype(t), (16, 16, 4, 16, 16, 4))

    # T1[ixu, ixd, jxu, jxd, iyu, iyd, jyu, jyd, iτu, iτd, ixup, ixdp, jxup, jxdp, iyup, iydp, jyup, jydp, iτup, iτdp]
    for ixu in 0:1, ixd in 0:1, jxu in 0:1, jxd in 0:1, iyu in 0:1, iyd in 0:1, jyu in 0:1, jyd in 0:1, iτu in 0:1, iτd in 0:1
        for ixup in 0:1, ixdp in 0:1, jxup in 0:1, jxdp in 0:1, iyup in 0:1, iydp in 0:1, jyup in 0:1, jydp in 0:1, iτup in 0:1, iτdp in 0:1

            Ix = f(ixu, ixd, jxu, jxd); Iy = f(iyu, iyd, jyu, jyd); Iτ = f(iτu, iτd)
            Ixp = f(ixup, ixdp, jxup, jxdp); Iyp = f(iyup, iydp, jyup, jydp); Iτp = f(iτup, iτdp)

            exponent1 = ixu + ixd + iyu + iyd
            exponent2 = ixu + jxu + ixup + jxup + ixd + jxd + ixdp + jxdp + iyu + jyu + iyup + jyup + iyd + jyd + iydp + jydp 
            exponent3 = ixu * iτu + ixd * (iτd + jyu + iyup + jxu + ixup + iτup + jyup + iyu + jxup + iτu) + jxu * (ixup + iτup + jyup + iyu + jxup + iτu) 
            + jxd * (ixdp + iτdp + jydp + iyd + jxdp + iτd + jyu + iyup + ixup + iτup + jyup + iyu + jxup + iτu) + iyu * (jxup + iτu)
            + iyd * (jxdp + iτd + jyu + iyup + ixup + iτup + jyup + jxup + iτu) + jyu * (iyup + ixup + iτup + jyup + jxup + iτu) + jyd * (iydp + ixdp + iτdp + 
            jydp + jxdp + iτd + iyup + ixup + iτup + jyup + jxup + iτu) + iτd * (iyup + ixup + iτup + jyup + jxup) + iτdp * (jydp + jxdp + iyup + ixup + iτup + jyup + jxup)
            + iτup * (jyup + jxup) + jydp * (jxdp + iyup + ixup + jyup + jxup) + jyup * jxup + iydp * (ixdp + jxdp + iyup + ixup + jxup) + iyup * (ixup + jxup) + jxdp * (ixup + jxup) 
            + ixdp * ixup
            sign_num = (-1)^exponent3
            # sign_num = auto_sign((iτu, ixu, jxup, iyu, jyup, iτup, ixup, jxu, iyup, jyu, iτd, ixd, jxdp, iyd, jydp, iτdp, ixdp, jxd, iydp, jyd), (2, 12, 8, 18, 4, 14, 10, 20, 1, 11, 16, 6, 15, 5, 19, 9, 13, 3, 17, 7))

            T_coef[Ix, Iy, Iτ, Ixp, Iyp, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * (
                isequal(iτd + ixd + jxdp + iyd + jydp, 1) * isequal(iτdp + ixdp + jxd + iydp + jyd, 1) * isequal(iτu + ixu + jxup + iyu + jyup, 1) * isequal(iτup + ixup + jxu + iyup + jyu, 1)
                - (μ*ϵ+1) * isequal(iτd + ixd + jxdp + iyd + jydp, 0) * isequal(iτdp + ixdp + jxd + iydp + jyd, 0) * isequal(iτu + ixu + jxup + iyu + jyup, 1) * isequal(iτup + ixup + jxu + iyup + jyu, 1)
                - (μ*ϵ+1) * isequal(iτd + ixd + jxdp + iyd + jydp, 1) * isequal(iτdp + ixdp + jxd + iydp + jyd, 1) * isequal(iτu + ixu + jxup + iyu + jyup, 0) * isequal(iτup + ixup + jxu + iyup + jyu, 0)
                - (U*ϵ-(μ*ϵ+1)^2) * isequal(iτd + ixd + jxdp + iyd + jydp, 0) * isequal(iτdp + ixdp + jxd + iydp + jyd, 0) * isequal(iτu + ixu + jxup + iyu + jyup, 0) * isequal(iτup + ixup + jxu + iyup + jyu, 0)
            )
             
        end
    end

    return T_coef
end

function Hubbard_model_Honeycomb_A(t, μ, ϵ, U)

    # T_coef[(iτu, iτd), (i1u, i1d, j1u, j1d), (i2u, i2d, j2u, j2d), (i3u, i3d, j3u, j3d), (iτup, iτdp)]
    T_coef = zeros(eltype(t), (4, 16, 16, 16, 4))

    for iτu in 0:1, iτd in 0:1, i1u in 0:1, i1d in 0:1, j1u in 0:1, j1d in 0:1, i2u in 0:1, i2d in 0:1, j2u in 0:1, j2d in 0:1
        for i3u in 0:1, i3d in 0:1, j3u in 0:1, j3d in 0:1, iτup in 0:1, iτdp in 0:1

            Iτ = f(iτu, iτd); Iτp = f(iτup, iτdp)
            I1 = f(i1u, i1d, j1u, j1d)
            I2 = f(i2u, i2d, j2u, j2d)
            I3 = f(i3u, i3d, j3u, j3d)
            
            exponent1 = iτu + iτd + iτup + iτdp + j1u + j1d + j2u + j2d + j3u + j3d
            exponent2 = i1u + i1d + j1u + j1d + i2u + i2d + j2u + j2d + i3u + i3d + j3u + j3d
            exponent3 = iτu * (iτu + i3d + i2d + i1d + iτd) + iτd * (iτd + iτu) + i1u * (i1u + i3d + i2d + i1d + iτu + iτd) +
            i1d * (i1d + i1u + iτu + iτd) +  j1u * (j2u + j3u) + j1d * (j2d + j3d + iτup + j2u + j3u) + 
            i2u * (i3u + iτdp + j2d + j3d + iτup + j2u + j3u) + i2d * (i3d + i3u + iτdp + j2d + j3d + iτup + j2u + j3u) + 
            j2u * j3u + j2d * (j3d + iτup + j3u) + i3u * (iτdp + j3d + iτup + j3u) + i3d * (iτdp + j3d + iτup + j3u) + j3d * iτup + iτdp * iτup
            sign_num = (-1)^exponent3
            # sign_num = auto_sign((j3u, j2u, j1u, iτup, j3d, j2d, j1d, iτdp, i3u, i2u, i1u, iτu, i3d, i2d, i1d, iτd), (12, 16, 11, 15, 3, 7, 10, 14, 2, 6, 9, 13, 1, 5, 8, 4))

            T_coef[Iτ, I1, I2, I3, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * (
                isequal(iτd + i1d + i2d + i3d, 1) * isequal(iτu + i1u + i2u + i3u, 1) * isequal(iτdp + j1d + j2d + j3d, 1) * isequal(iτup + j1u + j2u + j3u, 1)
                - (μ*ϵ+1) * isequal(iτd + i1d + i2d + i3d, 1) * isequal(iτu + i1u + i2u + i3u, 0) * isequal(iτdp + j1d + j2d + j3d, 1) * isequal(iτup + j1u + j2u + j3u, 0)
                - (μ*ϵ+1) * isequal(iτd + i1d + i2d + i3d, 0) * isequal(iτu + i1u + i2u + i3u, 1) * isequal(iτdp + j1d + j2d + j3d, 0) * isequal(iτup + j1u + j2u + j3u, 1)
                + (U*ϵ-(μ*ϵ+1)^2) * isequal(iτd + i1d + i2d + i3d, 0) * isequal(iτu + i1u + i2u + i3u, 0) * isequal(iτdp + j1d + j2d + j3d, 0) * isequal(iτup + j1u + j2u + j3u, 0)
            )
            
        end
    end

    return T_coef
end

function Hubbard_model_Honeycomb_B(t, μ, ϵ, U)

    # T_coef[(iτu, iτd), (i3up, i3dp, j3up, j3dp), (i2up, i2dp, j2up, j2dp), (i1up, i1dp, j1up, j1dp), (iτup, iτdp)]
    T_coef = zeros(eltype(t), (4, 16, 16, 16, 4))

    for iτu in 0:1, iτd in 0:1, i3up in 0:1, i3dp in 0:1, j3up in 0:1, j3dp in 0:1, i2up in 0:1, i2dp in 0:1, j2up in 0:1, j2dp in 0:1
        for i1up in 0:1, i1dp in 0:1, j1up in 0:1, j1dp in 0:1, iτup in 0:1, iτdp in 0:1

            Iτ = f(iτu, iτd); Iτp = f(iτup, iτdp)
            I1p = f(i1up, i1dp, j1up, j1dp)
            I2p = f(i2up, i2dp, j2up, j2dp)
            I3p = f(i3up, i3dp, j3up, j3dp)

            exponent1 = iτu + iτd + j1dp + j1up + j2dp + j2up + j3dp + j3up
            exponent2 = i1up + i1dp + j1up + j1dp + i2up + i2dp + j2up + j2dp + i3up + i3dp + j3up + j3dp
            exponent3 = iτu * (iτu + j3dp + j2dp + j1dp + iτd) + iτd * (iτd + iτu) + j3dp * (j3dp + j2dp + j1dp + iτd + iτu) + 
            j3up * (j3up + j2up + j1up + j2dp + j1dp + iτd + iτu + j3dp) + i3dp * (iτup + i1up + i2up + i3up) + 
            j2dp * (j2dp + j1dp + iτd + iτu + j3dp + j3up + i3dp + i3up) + j2up * (iτdp + i1dp + i2dp + iτup + i1up + i2up) + 
            i2dp * (iτup + i1up + i2up) + j1dp * (j1up + iτdp + i1dp + iτup + i1up) + j1up * (iτdp + i1dp + iτup + i1up) +
            i1dp * (iτup + i1up) + iτdp * iτup
            sign_num = (-1)^exponent3
            # sign_num = auto_sign((i3up, i2up, i1up, iτup, i3dp, i2dp, i1dp, iτdp, j3up, j2up, j1up, iτu, j3dp, j2dp, j1dp, iτd), (12, 16, 13, 9, 5, 1, 14, 10, 6, 2, 15, 11, 7, 3, 8, 4))

            T_coef[Iτ, I3p, I2p, I1p, Iτp] = (-1)^exponent1 * (sqrt(ϵ * t))^exponent2 * sign_num * (
                isequal(iτd + j1dp + j2dp + j3dp, 1) * isequal(iτu + j1up + j2up + j3up, 1) * isequal(iτdp + i1dp + i2dp + i3dp, 1) * isequal(iτup + i1up + i2up + i3up, 1)
                - (μ*ϵ+1) * isequal(iτd + j1dp + j2dp + j3dp, 1) * isequal(iτu + j1up + j2up + j3up, 0) * isequal(iτdp + i1dp + i2dp + i3dp, 1) * isequal(iτup + i1up + i2up + i3up, 0)
                - (μ*ϵ+1) * isequal(iτd + j1dp + j2dp + j3dp, 0) * isequal(iτu + j1up + j2up + j3up, 1) * isequal(iτdp + i1dp + i2dp + i3dp, 0) * isequal(iτup + i1up + i2up + i3up, 1)
                + (U*ϵ-(μ*ϵ+1)^2) * isequal(iτd + j1dp + j2dp + j3dp, 0) * isequal(iτu + j1up + j2up + j3up, 0) * isequal(iτdp + i1dp + i2dp + i3dp, 0) * isequal(iτup + i1up + i2up + i3up, 0)
            )
        end
    end

    return T_coef
end

########################### Gross-Neveu model on the square lattice ###########################

"""
Exact solutions:

    lnZ/V = 1.4515448845652446  at m=0
    lnZ/V = 0.84112502060056715 at m=-1
"""

function Gross_Neveu_model(μ, m, g2; r=1)

    # Ā[i1, i2, j1p, j2p]
    Ā = zeros(ComplexF64, (2, 2, 2, 2))

    for i1 in 0:1, i2 in 0:1, j1p in 0:1, j2p in 0:1

        if (i1, i2, j1p, j2p) == (1, 1, 0, 0)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - 1 + im
        elseif (i1, i2, j1p, j2p) == (1, 0, 1, 0)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - 1 - 1
        elseif (i1, i2, j1p, j2p) == (1, 0, 0, 1)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - 1 - im
        elseif (i1, i2, j1p, j2p) == (0, 1, 1, 0)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - im - 1
        elseif (i1, i2, j1p, j2p) == (0, 1, 0, 1)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - im - im
        elseif (i1, i2, j1p, j2p) == (0, 0, 1, 1)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = 1 - im
        end

    end

    # A[j1, j2, i1p, i2p]
    A = zeros(ComplexF64, (2, 2, 2, 2))

    for j1 in 0:1, j2 in 0:1, i1p in 0:1, i2p in 0:1

        if (j1, j2, i1p, i2p) == (1, 1, 0, 0)
            A[j1+1, j2+1, i1p+1, i2p+1] = 1 + im
        elseif (j1, j2, i1p, i2p) == (1, 0, 1, 0)
            A[j1+1, j2+1, i1p+1, i2p+1] = 1 + 1
        elseif (j1, j2, i1p, i2p) == (1, 0, 0, 1)
            A[j1+1, j2+1, i1p+1, i2p+1] = 1 - im
        elseif (j1, j2, i1p, i2p) == (0, 1, 1, 0)
            A[j1+1, j2+1, i1p+1, i2p+1] = - im + 1
        elseif (j1, j2, i1p, i2p) == (0, 1, 0, 1)
            A[j1+1, j2+1, i1p+1, i2p+1] = - im - im
        elseif (j1, j2, i1p, i2p) == (0, 0, 1, 1)
            A[j1+1, j2+1, i1p+1, i2p+1] = - 1 - im
        end
    end

    T_coef = zeros(ComplexF64, (4, 4, 4, 4))

    for i1 in 0:1, j1 in 0:1, i2 in 0:1, j2 in 0:1
        for i1p in 0:1, j1p in 0:1, i2p in 0:1, j2p in 0:1

            I = f(i1, j1); J = f(i2, j2)
            Ip = f(i1p, j1p); Jp = f(i2p, j2p)

            sign1 = i1 * (j1 + j2 + i1p + i2p) + i2 * (j2 + i1p + i2p) + j1p * (i1p + i2p) + j2p * i2p + i1p + i2p
            sign2 = i2 - j2 + i2p -j2p
            sign3 = i1 + j1 + i2 + j2 + i1p + j1p + i2p + j2p
            sign4 = i1 + i2 + j2 + i1p
            sign5 = i2 + j2 + i2p + j2p

            T_coef[I, J, Ip, Jp] = (-1)^sign1 * exp(0.5*μ*sign2) * (1/sqrt(2))^sign3 * (
                ((m+2*r)^2 + 2*g2) * isequal(i1 + i2 + j1p + j2p, 0) * isequal(j1 + j2 + i1p + i2p, 0) -
                (m+2*r) * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1) - 
                (-1)^sign4 * (+im)^sign5 * (m+2*r) * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1) - 
                Ā[i1+1, i2+1, j1p+1, j2p+1] * A[j1+1, j2+1, i1p+1, i2p+1]
            )

        end
    end

    return T_coef
end

function Gross_Neveu_model_chiral_condensate(μ, m; r=1)

    T_coef = zeros(ComplexF64, (4, 4, 4, 4))

    for i1 in 0:1, j1 in 0:1, i2 in 0:1, j2 in 0:1
        for i1p in 0:1, j1p in 0:1, i2p in 0:1, j2p in 0:1

            I = f(i1, j1); J = f(i2, j2)
            Ip = f(i1p, j1p); Jp = f(i2p, j2p)

            sign1 = i1 * (j1 + j2 + i1p + i2p) + i2 * (j2 + i1p + i2p) + j1p * (i1p + i2p) + j2p * i2p + i1p + i2p
            sign2 = i2 - j2 + i2p -j2p
            sign3 = i1 + j1 + i2 + j2 + i1p + j1p + i2p + j2p
            sign4 = i1 + i2 + j2 + i1p
            sign5 = i2 + j2 + i2p + j2p

            T_coef[I, J, Ip, Jp] = (-1)^sign1 * exp(0.5*μ*sign2) * (1/sqrt(2))^sign3 * (
                (-1)^sign4 * (+im)^sign5 * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1)
                + isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1)
                - 2 * (m+2*r) * isequal(i1 + i2 + j1p + j2p, 0) * isequal(j1 + j2 + i1p + i2p, 0)
                )
        end
    end

    return T_coef
end

function Gross_Neveu_model_pseudoscalar_condensate(μ; r=1)

    T_coef = zeros(ComplexF64, (4, 4, 4, 4))

    for i1 in 0:1, j1 in 0:1, i2 in 0:1, j2 in 0:1
        for i1p in 0:1, j1p in 0:1, i2p in 0:1, j2p in 0:1

            I = f(i1, j1); J = f(i2, j2)
            Ip = f(i1p, j1p); Jp = f(i2p, j2p)

            sign1 = i1 * (j1 + j2 + i1p + i2p) + i2 * (j2 + i1p + i2p) + j1p * (i1p + i2p) + j2p * i2p + i1p + i2p
            sign2 = i2 - j2 + i2p -j2p
            sign3 = i1 + j1 + i2 + j2 + i1p + j1p + i2p + j2p
            sign4 = i1 + i2 + j2 + i1p
            sign5 = i2 + j2 + i2p + j2p

            T_coef[I, J, Ip, Jp] = (-1)^sign1 * exp(0.5*μ*sign2) * (1/sqrt(2))^sign3 * im * (
                (-1)^sign4 * (+im)^sign5 * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1)
                - isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1))
        end
    end

    return T_coef
end

function Gross_Neveu_model_pseudoscalar_condensate_square(μ; r=1)

    T_coef = zeros(ComplexF64, (4, 4, 4, 4))

    for i1 in 0:1, j1 in 0:1, i2 in 0:1, j2 in 0:1
        for i1p in 0:1, j1p in 0:1, i2p in 0:1, j2p in 0:1

            I = f(i1, j1); J = f(i2, j2)
            Ip = f(i1p, j1p); Jp = f(i2p, j2p)

            sign1 = i1 * (j1 + j2 + i1p + i2p) + i2 * (j2 + i1p + i2p) + j1p * (i1p + i2p) + j2p * i2p + i1p + i2p
            sign2 = i2 - j2 + i2p -j2p
            sign3 = i1 + j1 + i2 + j2 + i1p + j1p + i2p + j2p
            sign4 = i1 + i2 + j2 + i1p
            sign5 = i2 + j2 + i2p + j2p

            T_coef[I, J, Ip, Jp] = (-1)^sign1 * exp(0.5*μ*sign2) * (1/sqrt(2))^sign3  * 2 *
            isequal(i1 + i2 + j1p + j2p, 0) * isequal(j1 + j2 + i1p + i2p, 0)
        end
    end

    return T_coef
end

function Gross_Neveu_model_external_field(μ, m, g2, h; r=1)

    # Ā[i1, i2, j1p, j2p]
    Ā = zeros(ComplexF64, (2, 2, 2, 2))

    for i1 in 0:1, i2 in 0:1, j1p in 0:1, j2p in 0:1

        if (i1, i2, j1p, j2p) == (1, 1, 0, 0)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - 1 + im
        elseif (i1, i2, j1p, j2p) == (1, 0, 1, 0)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - 1 - 1
        elseif (i1, i2, j1p, j2p) == (1, 0, 0, 1)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - 1 - im
        elseif (i1, i2, j1p, j2p) == (0, 1, 1, 0)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - im - 1
        elseif (i1, i2, j1p, j2p) == (0, 1, 0, 1)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = - im - im
        elseif (i1, i2, j1p, j2p) == (0, 0, 1, 1)
            Ā[i1+1, i2+1, j1p+1, j2p+1] = 1 - im
        end

    end

    # A[j1, j2, i1p, i2p]
    A = zeros(ComplexF64, (2, 2, 2, 2))

    for j1 in 0:1, j2 in 0:1, i1p in 0:1, i2p in 0:1

        if (j1, j2, i1p, i2p) == (1, 1, 0, 0)
            A[j1+1, j2+1, i1p+1, i2p+1] = 1 + im
        elseif (j1, j2, i1p, i2p) == (1, 0, 1, 0)
            A[j1+1, j2+1, i1p+1, i2p+1] = 1 + 1
        elseif (j1, j2, i1p, i2p) == (1, 0, 0, 1)
            A[j1+1, j2+1, i1p+1, i2p+1] = 1 - im
        elseif (j1, j2, i1p, i2p) == (0, 1, 1, 0)
            A[j1+1, j2+1, i1p+1, i2p+1] = - im + 1
        elseif (j1, j2, i1p, i2p) == (0, 1, 0, 1)
            A[j1+1, j2+1, i1p+1, i2p+1] = - im - im
        elseif (j1, j2, i1p, i2p) == (0, 0, 1, 1)
            A[j1+1, j2+1, i1p+1, i2p+1] = - 1 - im
        end
    end

    T_coef = zeros(ComplexF64, (4, 4, 4, 4))

    for i1 in 0:1, j1 in 0:1, i2 in 0:1, j2 in 0:1
        for i1p in 0:1, j1p in 0:1, i2p in 0:1, j2p in 0:1

            I = f(i1, j1); J = f(i2, j2)
            Ip = f(i1p, j1p); Jp = f(i2p, j2p)

            sign1 = i1 * (j1 + j2 + i1p + i2p) + i2 * (j2 + i1p + i2p) + j1p * (i1p + i2p) + j2p * i2p + i1p + i2p
            sign2 = i2 - j2 + i2p -j2p
            sign3 = i1 + j1 + i2 + j2 + i1p + j1p + i2p + j2p
            sign4 = i1 + i2 + j2 + i1p
            sign5 = i2 + j2 + i2p + j2p

            T_coef[I, J, Ip, Jp] = (-1)^sign1 * exp(0.5*μ*sign2) * (1/sqrt(2))^sign3 * (
                ((m+2*r)^2 + 2*g2 + h^2) * isequal(i1 + i2 + j1p + j2p, 0) * isequal(j1 + j2 + i1p + i2p, 0) -
                (m+2*r+im*h) * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1) - 
                (-1)^sign4 * (+im)^sign5 * (m+2*r-im*h) * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1) - 
                Ā[i1+1, i2+1, j1p+1, j2p+1] * A[j1+1, j2+1, i1p+1, i2p+1]
            )

        end
    end

    return T_coef
end

function Gross_Neveu_model_external_field_pseudoscalar(μ, h; r=1)

    T_coef = zeros(ComplexF64, (4, 4, 4, 4))

    for i1 in 0:1, j1 in 0:1, i2 in 0:1, j2 in 0:1
        for i1p in 0:1, j1p in 0:1, i2p in 0:1, j2p in 0:1

            I = f(i1, j1); J = f(i2, j2)
            Ip = f(i1p, j1p); Jp = f(i2p, j2p)

            sign1 = i1 * (j1 + j2 + i1p + i2p) + i2 * (j2 + i1p + i2p) + j1p * (i1p + i2p) + j2p * i2p + i1p + i2p
            sign2 = i2 - j2 + i2p -j2p
            sign3 = i1 + j1 + i2 + j2 + i1p + j1p + i2p + j2p
            sign4 = i1 + i2 + j2 + i1p
            sign5 = i2 + j2 + i2p + j2p

            T_coef[I, J, Ip, Jp] = (-1)^sign1 * exp(0.5*μ*sign2) * (1/sqrt(2))^sign3 * (
                im * (-1)^sign4 * (+im)^sign5 * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1)
                - im * isequal(i1 + i2 + j1p + j2p, 1) * isequal(j1 + j2 + i1p + i2p, 1)
                + 2 * h * isequal(i1 + i2 + j1p + j2p, 0) * isequal(j1 + j2 + i1p + i2p, 0)
                )

        end
    end

    return T_coef
end

########################### Two-flavor Gross-Neveu-Wilson model on the square lattice ###########################

"""
function Gross_Neveu_Wilson_model(μ, m, gσ2, gπ2; r=1)

    T_coef = zeros(ComplexF64, (16, 16, 16, 16))
    
    for i11 in 0:1, i12 in 0:1, j11 in 0:1, j12 in 0:1
        for i21 in 0:1, i22 in 0:1, j21 in 0:1, j22 in 0:1
            for i11p in 0:1, i12p in 0:1, j11p in 0:1, j12p in 0:1
                for i21p in 0:1, i22p in 0:1, j21p in 0:1, j22p in 0:1

                    I = f(i11, i12, j11, j12)
                    J = f(i21, i22, j21, j22)
                    Ip = f(i11p, i12p, j11p, j12p)
                    Jp = f(i21p, i22p, j21p, j22p)

                    sign1 = i11 + i12 + j11 + j12 + i21 + i22 + j21 + j22 + i11p + i12p + j11p + j12p + i21p + i22p + j21p + j22p
                    sign2 = i22p + i21p + i12p + i11p
                    sign3 = i21 + i21p + i22 + i22p - j21p - j21 - j22p - j22
                    sign_number = auto_sign(
                        (i11, i12, j11, j12, i21, i22, j21, j22, i11p, i12p, j11p, j12p, i21p, i22p, j21p, j22p), 
                        (1, 5, 9, 13, 2, 6, 10, 14, 7, 3, 15, 11, 8, 4, 16, 12)
                        )
                    sign4a = (i11p + j11) * (j11p + i11) + i11p + i21p + i11 + j21p
                    sign4b = i21p + j21 + j21p + i21
                    sign5a = (i12p + j12) * (j22p + i12) + i12p + i22p + i12 + j22p
                    sign5b = i12p + j12 + j12p + i12

                    coef1 = (m + 2*r)^2 * (3*gσ2 + gπ2) + (m + 2*r)^4 + gσ2 * gπ2/2
                    coef2 = A(j11, j21, i11p, i21p) * Ā(i11, i21, j11p, j21p) * ((m + 2*r)^2 + gσ2/2 + gπ2/2)
                    coef3 = A(j12, j22, i12p, i22p) * Ā(i12, i22, j12p, j22p) * ((m + 2*r)^2 + gσ2/2 + gπ2/2)
                    coef4 = (m + 2*r) * (3 * gσ2/2 + gπ2/2 + (m + 2*r)^2) * (1 - (-1)^sign4a * (im)^sign4b)
                    coef5 = (m + 2*r) * (3 * gσ2/2 + gπ2/2 + (m + 2*r)^2) * (1 - (-1)^sign5a * (im)^sign5b)
                    coef6 = (1 + (-1)^sign4a * (im)^sign4b * (-1)^sign5a * (im)^sign5b) * (gσ2/2 - gπ2/2) + 
                    ((-1)^sign4a * (im)^sign4b + (-1)^sign5a * (im)^sign5b) * (gσ2/2 + gπ2/2)
                    coef7 = A(j11, j21, i11p, i21p) * Ā(i11, i21, j11p, j21p) * A(j12, j22, i12p, i22p) * Ā(i12, i22, j12p, j22p)

                    T_coef[I, J, Ip, Jp] = 
                    coef1 * isequal(i11 + i21 + j11p + j21p, 0) * isequal(i12 + i22 + j12p + j22p, 0) * isequal(j11 + j21 + i11p + i21p, 0) * isequal(j12 + j22 + i12p + i22p, 0) -
                    coef2 * isequal(i11 + i21 + j11p + j21p, 2) * isequal(i12 + i22 + j12p + j22p, 0) * isequal(j11 + j21 + i11p + i21p, 2) * isequal(j12 + j22 + i12p + i22p, 0) - 
                    coef3 * isequal(i11 + i21 + j11p + j21p, 1) * isequal(i12 + i22 + j12p + j22p, 0) * isequal(j11 + j21 + i11p + i21p, 1) * isequal(j12 + j22 + i12p + i22p, 0) +
                    coef4 * isequal(i11 + i21 + j11p + j21p, 1) * isequal(i12 + i22 + j12p + j22p, 0) * isequal(j11 + j21 + i11p + i21p, 1) * isequal(j12 + j22 + i12p + i22p, 0) +
                    coef5 * isequal(i11 + i21 + j11p + j21p, 0) * isequal(i12 + i22 + j12p + j22p, 1) * isequal(j11 + j21 + i11p + i21p, 0) * isequal(j12 + j22 + i12p + i22p, 1) +
                    coef6 * isequal(i11 + i21 + j11p + j21p, 1) * isequal(i12 + i22 + j12p + j22p, 1) * isequal(j11 + j21 + i11p + i21p, 1) * isequal(j12 + j22 + i12p + i22p, 1) +
                    coef7 * isequal(i11 + i21 + j11p + j21p, 2) * isequal(i12 + i22 + j12p + j22p, 2) * isequal(j11 + j21 + i11p + i21p, 2) * isequal(j12 + j22 + i12p + i22p, 2)
                end
            end
        end
    end

end
"""

function Ā(i1, i2, j1p, j2p)
    if (i1, i2, j1p, j2p) == (1, 1, 0, 0)
        -1 + im
    elseif (i1, i2, j1p, j2p) == (1, 0, 1, 0)
        - 1 - 1
    elseif (i1, i2, j1p, j2p) == (1, 0, 0, 1)
        - 1 - im
    elseif (i1, i2, j1p, j2p) == (0, 1, 1, 0)
        - im - 1
    elseif (i1, i2, j1p, j2p) == (0, 1, 0, 1)
        - im - im
    elseif (i1, i2, j1p, j2p) == (0, 0, 1, 1)
        1 - im
    else
        0 + 0im
    end
end

function A(j1, j2, i1p, i2p)
    if (j1, j2, i1p, i2p) == (1, 1, 0, 0)
        1 + im
    elseif (j1, j2, i1p, i2p) == (1, 0, 1, 0)
        1 + 1
    elseif (j1, j2, i1p, i2p) == (1, 0, 0, 1)
        1 - im
    elseif (j1, j2, i1p, i2p) == (0, 1, 1, 0)
        - im + 1
    elseif (j1, j2, i1p, i2p) == (0, 1, 0, 1)
        - im - im
    elseif (j1, j2, i1p, i2p) == (0, 0, 1, 1)
        - 1 - im
    else
        0 + 0im
    end
end

P(a, b, c, d) = (-1)^(a + b + c + d)

I(a, b, c, d) = (im)^(a + b + c + d)

function Gross_Neveu_Wilson_model_Nf2(μ1, μ2, m1, m2, gσ2, gπ2, H)

    # T_coef[(i1, j1, k1, l1), (i2, j2, k2, l2), (i1p, j1p, k1p, l1p), (i2p, j2p, k2p, l2p)]
    T_coef = zeros(ComplexF64, (16, 16, 16, 16))

    C12 = -(m1 + 2 + im * H)
    C22 = -(m1 + 2 - im * H)
    C32 = -(m2 + 2 - im * H)
    C42 = -(m2 + 2 + im * H)

    C14 = (m1 + 2 + im * H) * (m2 + 2 - im * H) + (gσ2 - gπ2) / 2
    C24 = (m1 + 2 + im * H) * (m2 + 2 + im * H) + (gσ2 + gπ2) / 2
    C34 = (m2 + 2 - im * H) * (m1 + 2 - im * H) + (gσ2 + gπ2) / 2
    C44 = (m2 + 2 + im * H) * (m1 + 2 - im * H) + (gσ2 - gπ2) / 2
    C54 = (m1 + 2 + im * H) * (m1 + 2 - im * H) + (gσ2 + gπ2) / 2
    C64 = (m2 + 2 - im * H) * (m2 + 2 + im * H) + (gσ2 + gπ2) / 2

    C16 = -(m1 + 2 + im * H) * (m1 + 2 - im * H) * (m2 + 2 - im * H) -
          (gσ2 + gπ2)/2 * (m1 + 2 + im * H) - 
          (gσ2 - gπ2)/2 * (m1 + 2 - im * H) -
          (gσ2 + gπ2)/2 * (m2 + 2 - im * H)

    C26 = -(m1 + 2 + im * H) * (m1 + 2 - im * H) * (m2 + 2 + im * H) - 
          (gσ2 - gπ2)/2 * (m1 + 2 + im * H) - 
          (gσ2 + gπ2)/2 * (m1 + 2 - im * H) - 
          (gσ2 + gπ2)/2 * (m2 + 2 + im * H)

    C36 = -(m1 + 2 + im * H) * (m2 + 2 - im * H) * (m2 + 2 + im * H) -
          (gσ2 + gπ2)/2 * (m1 + 2 + im * H) - 
          (gσ2 + gπ2)/2 * (m2 + 2 - im * H) - 
          (gσ2 - gπ2)/2 * (m2 + 2 + im * H)

    C46 = -(m1 + 2 - im * H) * (m2 + 2 - im * H) * (m2 + 2 + im * H) - 
          (gσ2 + gπ2)/2 * (m1 + 2 - im * H) - 
          (gσ2 - gπ2)/2 * (m2 + 2 - im * H) - 
          (gσ2 + gπ2)/2 * (m2 + 2 + im * H)

    C8 = (m1 + 2 + im * H) * (m1 + 2 - im * H) * (m2 + 2 + im * H) * (m2 + 2 - im * H) + 
         ((gσ2 - gπ2)/2)^2 + 2 * ((gσ2 + gπ2)/2)^2 + 
         (gσ2 - gπ2)/2 * (m1 + 2 - im * H) * (m2 + 2 + im * H) + 
         (gσ2 - gπ2)/2 * (m1 + 2 + im * H) * (m2 + 2 - im * H) + 
         (gσ2 + gπ2)/2 * (m2 + 2 - im * H) * (m2 + 2 + im * H) + 
         (gσ2 + gπ2)/2 * (m1 + 2 + im * H) * (m1 + 2 - im * H) + 
         (gσ2 + gπ2)/2 * (m1 + 2 - im * H) * (m2 + 2 - im * H) + 
         (gσ2 + gπ2)/2 * (m1 + 2 + im * H) * (m2 + 2 + im * H)

    for i1 in 0:1, j1 in 0:1, k1 in 0:1, l1 in 0:1
        for i2 in 0:1, j2 in 0:1, k2 in 0:1, l2 in 0:1
            for i1p in 0:1, j1p in 0:1, k1p in 0:1, l1p in 0:1
                for i2p in 0:1, j2p in 0:1, k2p in 0:1, l2p in 0:1

                    I1 = f(i1, j1, k1, l1)
                    I2 = f(i2, j2, k2, l2)
                    I1p = f(i1p, j1p, k1p, l1p)
                    I2p = f(i2p, j2p, k2p, l2p)

                    s1 = i1p + i2p + k1p + k2p
                    s2 = i2 - j2 + i2p - j2p
                    s3 = k2 - l2 + k2p - l2p
                    s4 = i1 + j1 + i2 + j2 + i1p + j1p + i2p + j2p
                    s5 = k1 + l1 + k2 + l2 + k1p + l1p + k2p + l2p

                    s6 = 
                    i1 * (l1 + l2 + k1p + k2p + k1 + k2 + l1p + l2p + j1 + j2 + i1p + i2p) + 
                    j1 * (l1 + l2 + k1p + k2p + k1 + k2 + l1p + l2p) + 
                    k1 * (l1 + l2 + k1p + k2p) + 
                    i2 * (l2 + k1p + k2p + k2 + l1p + l2p + j2 + i1p + i2p) + 
                    j2 * (l2 + k1p + k2p + k2 + l1p + l2p) + 
                    k2 * (l2 + k1p + k2p) + 
                    l1p * (k1p + k2p) + 
                    j1p * (k2p + l2p + i1p + i2p) + 
                    i1p * (k2p + l2p) + 
                    l2p * k2p + 
                    j2p * i2p

                    t0 = Ā(i1, i2, j1p, j2p) * A(j1, j2, i1p, i2p) * Ā(k1, k2, l1p, l2p) * A(l1, l2, k1p, k2p)

                    t2 = 
                    (-C12 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) - C22) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    Ā(k1, k2, l1p, l2p) *  
                    A(l1, l2, k1p, k2p) + 
                    (-C32 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) - C42) * 
                    Ā(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1)

                    t4 = 
                    (C14 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + 
                    C24 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) + 
                    C34 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + 
                    C44) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    C54 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    Ā(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p) - 
                    C64 * 
                    Ā(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    t6 = 
                    (C16 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + C26) * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) + 
                    (C36 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) + C46) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    t8 = C8 * isequal(j2p + j1p + i2 + i1, 0) * isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 0) * isequal(k2p + k1p + l2 + l1, 0)

                    t = t0 + t2 + t4 + t6 + t8

                    T_coef[I1, I2, I1p, I2p] = (-1)^(s1) * exp(0.5 * μ1 * s2) * exp(0.5 * μ2 * s3) * 
                    (1/sqrt(2))^s4 * (1/sqrt(2))^s5 * (-1)^s6 * t
                end
            end
        end
    end
    
    return T_coef
end

function Gross_Neveu_Wilson_model_Nf2_impurity(m1, m2, gσ2, gπ2, H)

    # T_coef[(i1, j1, k1, l1), (i2, j2, k2, l2), (i1p, j1p, k1p, l1p), (i2p, j2p, k2p, l2p)]
    T_coef_chiral_condensate = zeros(ComplexF64, (16, 16, 16, 16))
    T_coef_pseudoscalar_singlet = zeros(ComplexF64, (16, 16, 16, 16))
    T_coef_pseudoscalar_triplet = zeros(ComplexF64, (16, 16, 16, 16))

    a114 = -(m2 + 2 - im * H)   
    b114 = -(m2 + 2 + im * H)   
    c114 = -(m1 + 2 - im * H)

    a214 = -(m2 + 2 - im * H)  
    b214 = -(m2 + 2 + im * H)   
    c214 = -(m1 + 2 + im * H) 

    a124 = -(m1 + 2 + im * H)
    b124 = -(m1 + 2 - im * H)
    c124 = -(m2 + 2 + im * H) 

    a224 = -(m1 + 2 + im * H)
    b224 = -(m1 + 2 - im * H)   
    c224 = -(m2 + 2 - im * H)
    
    a116 = (m1 + 2 - im * H) * (m2 + 2 - im * H) + (gσ2 + gπ2)/2
    b116 = (m1 + 2 - im * H) * (m2 + 2 + im * H) + (gσ2 + gπ2)/2
    c116 = (m2 + 2 - im * H) * (m2 + 2 + im * H) + (gσ2 - gπ2)/2

    a216 = (m1 + 2 + im * H) * (m2 + 2 - im * H) + (gσ2 - gπ2)/2
    b216 = (m1 + 2 + im * H) * (m2 + 2 + im * H) + (gσ2 + gπ2)/2
    c216 = (m2 + 2 - im * H) * (m2 + 2 + im * H) + (gσ2 + gπ2)/2

    a126 = (m1 + 2 + im * H) * (m1 + 2 - im * H) + (gσ2 + gπ2)/2
    b126 = (m1 + 2 + im * H) * (m2 + 2 + im * H) + (gσ2 + gπ2)/2
    c126 = (m1 + 2 - im * H) * (m2 + 2 + im * H) + (gσ2 - gπ2)/2

    a226 = (m1 + 2 + im * H) * (m1 + 2 - im * H) + (gσ2 + gπ2)/2
    b226 = (m1 + 2 + im * H) * (m2 + 2 - im * H) + (gσ2 - gπ2)/2
    c226 = (m1 + 2 - im * H) * (m2 + 2 - im * H) + (gσ2 + gπ2)/2

    a118 = -(m1 + 2 - im * H) * (m2 + 2 - im * H) * (m2 + 2 + im * H) - 
    (m1 + 2 - im * H) * (gσ2 + gπ2)/2 - 
    (m2 + 2 - im * H) * (gσ2 - gπ2)/2 - 
    (m2 + 2 + im * H) * (gσ2 + gπ2)/2

    a218 = -(m1 + 2 + im * H) * (m2 + 2 - im * H) * (m2 + 2 + im * H) - 
    (m1 + 2 + im * H) * (gσ2 + gπ2)/2 - 
    (m2 + 2 - im * H) * (gσ2 + gπ2)/2 - 
    (m2 + 2 + im * H) * (gσ2 - gπ2)/2

    a128 = -(m1 + 2 + im * H) * (m1 + 2 - im * H) * (m2 + 2 + im * H) - 
    (m1 + 2 + im * H) * (gσ2 - gπ2)/2 - 
    (m1 + 2 - im * H) * (gσ2 + gπ2)/2 - 
    (m2 + 2 + im * H) * (gσ2 + gπ2)/2

    a228 = -(m1 + 2 + im * H) * (m1 + 2 - im * H) * (m2 + 2 - im * H) - 
    (m1 + 2 + im * H) * (gσ2 + gπ2)/2 - 
    (m1 + 2 - im * H) * (gσ2 - gπ2)/2 - 
    (m2 + 2 - im * H) * (gσ2 + gπ2)/2

    for i1 in 0:1, j1 in 0:1, k1 in 0:1, l1 in 0:1
        for i2 in 0:1, j2 in 0:1, k2 in 0:1, l2 in 0:1
            for i1p in 0:1, j1p in 0:1, k1p in 0:1, l1p in 0:1
                for i2p in 0:1, j2p in 0:1, k2p in 0:1, l2p in 0:1

                    I1 = f(i1, j1, k1, l1)
                    I2 = f(i2, j2, k2, l2)
                    I1p = f(i1p, j1p, k1p, l1p)
                    I2p = f(i2p, j2p, k2p, l2p)

                    s112 = -P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    Ā(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s212 = -isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    Ā(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s122 = -P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) * 
                    Ā(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1)

                    s222 = -Ā(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1)

                    s114 = (a114 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + 
                    b114 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2)) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    c114 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    Ā(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s214 = (a214 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + b214) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    c214 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    Ā(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s124 = (a124 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + 
                    b124 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2)) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    c124 * 
                    Ā(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0) 
                    
                    s224 = (a224 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) + b224) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    c224 * 
                    Ā(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s116 = (a116 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + b116) * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) + 
                    c116 * 
                    P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s216 = (a216 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + b216) * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) + 
                    c216 * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s126 = a126 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) + 
                    (b126 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) + c126) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s226 = a226 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) + 
                    (b226 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) + c226) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s118 = a118 * isequal(j2p + j1p + i2 + i1, 0) * isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 0) * isequal(k2p + k1p + l2 + l1, 0)

                    s218 = a218 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s128 = a128 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s228 = a228 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0)

                    s11 = s112 + s114 + s116 + s118
                    s21 = s212 + s214 + s216 + s218
                    s12 = s122 + s124 + s126 + s128
                    s22 = s222 + s224 + s226 + s228    

                    T_coef_chiral_condensate[I1, I2, I1p, I2p] = s11 + s21 + s12 + s22
                    T_coef_pseudoscalar_singlet[I1, I2, I1p, I2p] = im * s11 - im * s21 + im * s12 - im * s22
                    T_coef_pseudoscalar_triplet[I1, I2, I1p, I2p] = im * s11 - im * s21 - im * s12 + im * s22
                end
            end
        end
    end
    
    return T_coef_chiral_condensate, T_coef_pseudoscalar_singlet, T_coef_pseudoscalar_triplet
end

########################### Z2-symmetric Ising tensor for testing purposes ###########################

function isingtensor2d_forloop(β::Float64)
    
    W = [exp(β) exp(-β); exp(-β) exp(β)]
    U, S, _ = svd(W)
    Q = U * Diagonal(sqrt.(S))
    T = zeros(Float64, (2, 2, 2, 2))
    for i = 1:2, j = 1:2, k = 1:2, l = 1:2
        T[i, j, k, l] += Q[1, i] * Q[1, j] * Q[1, k] * Q[1, l]
        T[i, j, k, l] += Q[2, i] * Q[2, j] * Q[2, k] * Q[2, l]
    end
    
    return T
end

function isingtensor3d_forloop(β::Float64)
    
    W = [exp(β) exp(-β); exp(-β) exp(β)]
    U, S, _ = svd(W)
    Q = U * Diagonal(sqrt.(S))
    T = zeros(Float64, (2, 2, 2, 2, 2, 2))
    for i = 1:2, j = 1:2, k = 1:2, l = 1:2, m = 1:2, n = 1:2
        T[i, j, k, l, m, n] += Q[1, i] * Q[1, j] * Q[1, k] * Q[1, l] * Q[1, m] * Q[1, n]
        T[i, j, k, l, m, n] += Q[2, i] * Q[2, j] * Q[2, k] * Q[2, l] * Q[2, m] * Q[2, n]
    end
    
    return T
end 

########################### Basic (Z2-graded) index fusion functions ###########################

function f(ix::Int64)

    Ix = 0

    if ix == 0
        Ix = 1
    else
        Ix = 2
    end

    return Ix
end

function f(ix::Int64, jx::Int64)

    Ix = 0

    if (ix, jx) == (0, 0)
        Ix = 1
    elseif (ix, jx) == (1, 1)
        Ix = 2
    elseif (ix, jx) == (0, 1)
        Ix = 3
    else
        Ix = 4
    end

    return Ix
end

function f(ix::Int, jx::Int, kx::Int)

    Ix = 0

    if (ix, jx, kx) == (0, 0, 0)
        Ix = 1
    elseif (ix, jx, kx) == (1, 1, 0)
        Ix = 2
    elseif (ix, jx, kx) == (1, 0, 1)
        Ix = 3
    elseif (ix, jx, kx) == (0, 1, 1)
        Ix = 4
    elseif (ix, jx, kx) == (1, 0, 0)
        Ix = 5
    elseif (ix, jx, kx) == (0, 1, 0)
        Ix = 6
    elseif (ix, jx, kx) == (0, 0, 1)
        Ix = 7
    else
        Ix = 8
    end

    return Ix
end

function f(ixu::Int64, ixd::Int64, jxu::Int64, jxd::Int64)

    Ix = 0

    if (ixu, ixd, jxu, jxd) == (0, 0, 0, 0)
        Ix = 1
    elseif (ixu, ixd, jxu, jxd) == (1, 1, 0, 0)
        Ix = 2
    elseif (ixu, ixd, jxu, jxd) == (1, 0, 1, 0)
        Ix = 3
    elseif (ixu, ixd, jxu, jxd) == (1, 0, 0, 1)
        Ix = 4
    elseif (ixu, ixd, jxu, jxd) == (0, 1, 1, 0)
        Ix = 5
    elseif (ixu, ixd, jxu, jxd) == (0, 1, 0, 1)
        Ix = 6
    elseif (ixu, ixd, jxu, jxd) == (0, 0, 1, 1)
        Ix = 7
    elseif (ixu, ixd, jxu, jxd) == (1, 1, 1, 1)
        Ix = 8
    elseif (ixu, ixd, jxu, jxd) == (1, 0, 0, 0)
        Ix = 9
    elseif (ixu, ixd, jxu, jxd) == (0, 1, 0, 0)
        Ix = 10
    elseif (ixu, ixd, jxu, jxd) == (0, 0, 1, 0)
        Ix = 11
    elseif (ixu, ixd, jxu, jxd) == (1, 1, 1, 0)
        Ix = 12
    elseif (ixu, ixd, jxu, jxd) == (0, 0, 0, 1)
        Ix = 13
    elseif (ixu, ixd, jxu, jxd) == (1, 1, 0, 1)
        Ix = 14
    elseif (ixu, ixd, jxu, jxd) == (1, 0, 1, 1)
        Ix = 15
    else
        Ix = 16
    end

    return Ix
end

function f(ixu::Int64, ixd::Int64, jxu::Int64, jxd::Int64, kxu::Int64)

    Ix1 = f(ixu, ixd, jxu, jxd)
    Ix2 = f(kxu)

    if (Ix1 <= 8) && (Ix2 == 1)
        return Ix1
    elseif (Ix1 > 8) && (Ix2 == 2)
        return Ix1
    elseif (Ix1 <= 8) && (Ix2 == 2)
        return Ix1 + 16
    else
        return Ix1 + 16
    end
end
