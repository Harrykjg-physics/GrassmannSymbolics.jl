using Documenter
using GrassmannSymbolics

DocMeta.setdocmeta!(
    GrassmannSymbolics,
    :DocTestSetup,
    :(using GrassmannSymbolics);
    recursive=true,
)

makedocs(;
    modules=[GrassmannSymbolics],
    authors="GrassmannSymbolics Contributors",
    sitename="GrassmannSymbolics.jl",
    remotes=nothing,
    format=Documenter.HTML(;
        canonical="https://harrykjg-physics.github.io/GrassmannSymbolics/",
        repolink="https://github.com/Harrykjg-physics/GrassmannSymbolics",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Quick start" => "quickstart.md",
        "Grassmann algebra" => "algebra.md",
        "Local tensor compiler" => "local_tensors.md",
        "Coefficient arrays" => "coefficient_arrays.md",
        "Model examples" => "examples.md",
        "Registration notes" => "registration.md",
        "API index" => "api.md",
    ],
    checkdocs=:none,
)

if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo="github.com/Harrykjg-physics/GrassmannSymbolics.git",
        devbranch="main",
        push_preview=true,
    )
end
