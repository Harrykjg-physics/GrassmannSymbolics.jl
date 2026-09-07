using Pkg
const PROJECT_ROOT = dirname(@__DIR__)
Pkg.activate(PROJECT_ROOT)
using GrassmannSymbolics
include(joinpath(PROJECT_ROOT, "examples", "two_flavor_Gross_Neveu_Wilson.jl"))
include(get(ENV, "GRASSMANNTN_TENSOR_JL", joinpath(PROJECT_ROOT, "external_tensor.jl")))
include(joinpath(PROJECT_ROOT, "tmp", "new_nf2_impurity.jl"))
const F4 = [(0,0,0,0),(1,1,0,0),(1,0,1,0),(1,0,0,1),(0,1,1,0),(0,1,0,1),(0,0,1,1),(1,1,1,1),(1,0,0,0),(0,1,0,0),(0,0,1,0),(1,1,1,0),(0,0,0,1),(1,1,0,1),(1,0,1,1),(0,1,1,1)]
unpack(idx)=Tuple(vcat(collect(F4[idx[1]]),collect(F4[idx[2]]),collect(F4[idx[3]]),collect(F4[idx[4]])))
function correction(b)
    e=0
    for (a,c) in ((1,3),(1,4),(1,7),(1,8),(1,11),(1,12),(2,3),(2,4),(2,7),(2,8),(2,11),(2,12),(3,5),(3,6),(3,9),(3,10),(4,5),(4,6),(4,9),(4,10),(5,7),(5,8),(5,11),(5,12),(6,7),(6,8),(6,11),(6,12),(7,9),(7,10),(8,9),(8,10)) e ⊻= (b[a]&b[c]) end
    isodd(e) ? -1.0 : 1.0
end
function legs_direct(vars, idx)
    i1,j1,k1,l1,i2,j2,k2,l2,i1p,j1p,k1p,l1p,i2p,j2p,k2p,l2p=unpack(idx)
    ordered=GrassmannGenerator[vars.alpha[1],vars.beta[1],vars.alpha[2],vars.beta[2], vars.eta[1],vars.zeta[1],vars.eta[2],vars.zeta[2], vars.betabar[1],vars.alphabar[1],vars.betabar[2],vars.alphabar[2], vars.zetabar[1],vars.etabar[1],vars.zetabar[2],vars.etabar[2]]
    occ=(i1,j1,k1,l1, i2,j2,k2,l2, j1p,i1p,l1p,k1p, j2p,i2p,l2p,k2p)
    GrassmannGenerator[g for (g,b) in zip(ordered,occ) if b==1]
end
function numeric_model()
    vars=_tfgnw_generators(); params=(;M1=0.37,M2=-0.21,gs2=0.44,gp2=0.13,H=0.17,q=1/sqrt(2),z1=1.0,z2=1.0)
    base=_tfgnw_base_integrand(vars, params)
    n11=_tfgnw_density(vars,1,1); n12=_tfgnw_density(vars,1,2); n21=_tfgnw_density(vars,2,1); n22=_tfgnw_density(vars,2,2)
    s11=_tfgnw_integrate(vars,n11*base;simplify=false); s12=_tfgnw_integrate(vars,n12*base;simplify=false); s21=_tfgnw_integrate(vars,n21*base;simplify=false); s22=_tfgnw_integrate(vars,n22*base;simplify=false)
    (;vars,chiral=s11+s12+s21+s22,singlet=im*(s11-s12+s21-s22),triplet=im*(s11-s12-s21+s22))
end
function arr(tensor, vars)
    out=zeros(ComplexF64,16,16,16,16)
    for I in CartesianIndices(out)
        idx=Tuple(I); out[I]=correction(unpack(idx))*ComplexF64(get_coeff(tensor, legs_direct(vars,idx)))
    end
    out
end
model=numeric_model(); ch,ps,pt=Gross_Neveu_Wilson_model_Nf2_impurity(0.0,0.0,0.37,-0.21,0.44,0.13,0.17)
for (name, ours, ext) in (("chiral",arr(model.chiral,model.vars),ch),("singlet",arr(model.singlet,model.vars),ps),("triplet",arr(model.triplet,model.vars),pt))
    println("OBS ",name)
    n=0
    for I in CartesianIndices(ours)
        d=ours[I]-ext[I]
        if abs(d)>1e-10
            n+=1
            println("idx=",Tuple(I)," bits=",map(i->F4[i],Tuple(I))," diff=",d," diff_over_gp2half=",d/0.065," ours=",ours[I]," ext=",ext[I])
            n>=40 && break
        end
    end
    println("count_shown=",n," total=",count(abs.(ours-ext) .> 1e-10))
end
