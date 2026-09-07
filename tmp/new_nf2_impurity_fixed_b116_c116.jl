function Gross_Neveu_Wilson_model_Nf2_impurity(渭1, 渭2, m1, m2, g蟽2, g蟺2, H)

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
    
    a116 = (m1 + 2 - im * H) * (m2 + 2 - im * H) + (g蟽2 + g蟺2)/2
    b116 = (m1 + 2 - im * H) * (m2 + 2 + im * H) + (g蟽2 - g蟺2)/2
    c116 = (m2 + 2 - im * H) * (m2 + 2 + im * H) + (g蟽2 + g蟺2)/2

    a216 = (m1 + 2 + im * H) * (m2 + 2 - im * H) + (g蟽2 - g蟺2)/2
    b216 = (m1 + 2 + im * H) * (m2 + 2 + im * H) + (g蟽2 + g蟺2)/2
    c216 = (m2 + 2 - im * H) * (m2 + 2 + im * H) + (g蟽2 + g蟺2)/2

    a126 = (m1 + 2 + im * H) * (m1 + 2 - im * H) + (g蟽2 + g蟺2)/2
    b126 = (m1 + 2 + im * H) * (m2 + 2 + im * H) + (g蟽2 + g蟺2)/2
    c126 = (m1 + 2 - im * H) * (m2 + 2 + im * H) + (g蟽2 - g蟺2)/2

    a226 = (m1 + 2 + im * H) * (m1 + 2 - im * H) + (g蟽2 + g蟺2)/2
    b226 = (m1 + 2 + im * H) * (m2 + 2 - im * H) + (g蟽2 - g蟺2)/2
    c226 = (m1 + 2 - im * H) * (m2 + 2 - im * H) + (g蟽2 + g蟺2)/2

    a118 = -(m1 + 2 - im * H) * (m2 + 2 - im * H) * (m2 + 2 + im * H) - 
    (m1 + 2 - im * H) * (g蟽2 + g蟺2)/2 - 
    (m2 + 2 - im * H) * (g蟽2 - g蟺2)/2 - 
    (m2 + 2 + im * H) * (g蟽2 + g蟺2)/2

    a218 = -(m1 + 2 + im * H) * (m2 + 2 - im * H) * (m2 + 2 + im * H) - 
    (m1 + 2 + im * H) * (g蟽2 + g蟺2)/2 - 
    (m2 + 2 - im * H) * (g蟽2 + g蟺2)/2 - 
    (m2 + 2 + im * H) * (g蟽2 - g蟺2)/2

    a128 = -(m1 + 2 + im * H) * (m1 + 2 - im * H) * (m2 + 2 + im * H) - 
    (m1 + 2 + im * H) * (g蟽2 - g蟺2)/2 - 
    (m1 + 2 - im * H) * (g蟽2 + g蟺2)/2 - 
    (m2 + 2 + im * H) * (g蟽2 + g蟺2)/2

    a228 = -(m1 + 2 + im * H) * (m1 + 2 - im * H) * (m2 + 2 - im * H) - 
    (m1 + 2 + im * H) * (g蟽2 + g蟺2)/2 - 
    (m1 + 2 - im * H) * (g蟽2 - g蟺2)/2 - 
    (m2 + 2 - im * H) * (g蟽2 + g蟺2)/2

    for i1 in 0:1, j1 in 0:1, k1 in 0:1, l1 in 0:1
        for i2 in 0:1, j2 in 0:1, k2 in 0:1, l2 in 0:1
            for i1p in 0:1, j1p in 0:1, k1p in 0:1, l1p in 0:1
                for i2p in 0:1, j2p in 0:1, k2p in 0:1, l2p in 0:1

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

                    weight = (-1)^(s1) * exp(0.5 * 渭1 * s2) * exp(0.5 * 渭2 * s3) * 
                    (1/sqrt(2))^s4 * (1/sqrt(2))^s5 * (-1)^s6

                    I1 = f(i1, j1, k1, l1)
                    I2 = f(i2, j2, k2, l2)
                    I1p = f(i1p, j1p, k1p, l1p)
                    I2p = f(i2p, j2p, k2p, l2p)

                    s112 = -P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    膧(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s212 = -isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    膧(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s122 = -P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) * 
                    膧(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1)

                    s222 = -膧(i1, i2, j1p, j2p) * 
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
                    膧(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s214 = (a214 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + b214) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    c214 * 
                    isequal(j2p + j1p + i2 + i1, 0) * 
                    isequal(i2p + i1p + j2 + j1, 0) * 
                    膧(k1, k2, l1p, l2p) * 
                    A(l1, l2, k1p, k2p)

                    s124 = (a124 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2) + 
                    b124 * P(k2, k1, k1p, l2) * I(l2p, k2, k2p, l2)) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    c124 * 
                    膧(i1, i2, j1p, j2p) * 
                    A(j1, j2, i1p, i2p) * 
                    isequal(l2p + l1p + k2 + k1, 0) * 
                    isequal(k2p + k1p + l2 + l1, 0) 
                    
                    s224 = (a224 * P(i2, i1, i1p, j2) * I(j2p, i2, i2p, j2) + b224) * 
                    isequal(j2p + j1p + i2 + i1, 1) * 
                    isequal(i2p + i1p + j2 + j1, 1) * 
                    isequal(l2p + l1p + k2 + k1, 1) * 
                    isequal(k2p + k1p + l2 + l1, 1) - 
                    c224 * 
                    膧(i1, i2, j1p, j2p) * 
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

                    T_coef_chiral_condensate[I1, I2, I1p, I2p] = weight * (s11 + s21 + s12 + s22)
                    T_coef_pseudoscalar_singlet[I1, I2, I1p, I2p] = weight * (im * s11 - im * s21 + im * s12 - im * s22)
                    T_coef_pseudoscalar_triplet[I1, I2, I1p, I2p] = weight * (im * s11 - im * s21 - im * s12 + im * s22)
                end
            end
        end
    end
    
    return T_coef_chiral_condensate, T_coef_pseudoscalar_singlet, T_coef_pseudoscalar_triplet
end
