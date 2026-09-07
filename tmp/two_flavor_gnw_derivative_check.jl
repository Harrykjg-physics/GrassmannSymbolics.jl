using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
include(get(ENV, "GRASSMANNTN_TENSOR_JL", joinpath(PROJECT_ROOT, "external_tensor.jl")))
include(joinpath(PROJECT_ROOT, "tmp", "new_nf2_impurity.jl"))
function summarize(name, lhs, rhs)
    diff=lhs-rhs
    maxdiff=maximum(abs.(diff)); maxidx=argmax(abs.(diff)); ndiff=count(abs.(diff) .> 1e-7)
    println(name," maxdiff=",maxdiff," maxidx=",Tuple(maxidx)," lhs=",lhs[maxidx]," rhs=",rhs[maxidx]," ndiff=",ndiff)
end
function fd_chiral(mu1,mu2,m1,m2,gs2,gp2,H; h=1e-6)
    d1=(Gross_Neveu_Wilson_model_Nf2(mu1,mu2,m1+h,m2,gs2,gp2,H)-Gross_Neveu_Wilson_model_Nf2(mu1,mu2,m1-h,m2,gs2,gp2,H))/(2h)
    d2=(Gross_Neveu_Wilson_model_Nf2(mu1,mu2,m1,m2+h,gs2,gp2,H)-Gross_Neveu_Wilson_model_Nf2(mu1,mu2,m1,m2-h,gs2,gp2,H))/(2h)
    return -(d1+d2)
end
for gp2 in (0.13,0.0)
    println("GP2 ",gp2)
    fd=fd_chiral(0.0,0.0,0.37,-0.21,0.44,gp2,0.17)
    ch,ps,pt=Gross_Neveu_Wilson_model_Nf2_impurity(0.0,0.0,0.37,-0.21,0.44,gp2,0.17)
    summarize("external_new_chiral_vs_minus_mass_derivative", ch, fd)
end
