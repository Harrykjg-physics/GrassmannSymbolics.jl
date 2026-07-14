# Registration notes

This package is prepared as a standard Julia package:

```toml
name = "GrassmannSymbolics"
uuid = "4c7e4f2a-8b3d-4e9a-b5c1-3f2d1e0a9c7b"
version = "0.1.0"
```

## Standard Julia General workflow

The normal registration route is:

1. Ensure the package repository is public.
2. Ensure `Project.toml` has a valid `name`, `uuid`, `version`, `[deps]`, and
   `[compat]`.
3. Run tests locally or in CI.
4. Commit and push the release commit.
5. Trigger JuliaRegistrator on the release commit by commenting:

   ```text
   @JuliaRegistrator register
   ```

   on the commit, issue, or pull request.

6. Registrator opens a pull request to
   `https://github.com/JuliaRegistries/General`.
7. Wait for AutoMerge checks.  If AutoMerge does not merge automatically, respond
   to reviewer comments in the generated General registry pull request.
8. After the General registry pull request is merged, TagBot can create the
   corresponding GitHub tag and release.  Manual pre-registration tags are not
   required for the normal Registrator workflow.

Useful links:

- General registry: `https://github.com/JuliaRegistries/General`
- Registrator: `https://github.com/JuliaRegistries/Registrator.jl`
- Package naming guidelines:
  `https://pkgdocs.julialang.org/v1/creating-packages/`
- Package manager manual: `https://pkgdocs.julialang.org/v1/`

## Documentation deployment

This repository includes a Documenter.jl workflow.  Once GitHub Pages is enabled
for the generated `gh-pages` branch, the documentation URL is:

```text
https://harrykjg-physics.github.io/GrassmannSymbolics/
```

The workflow file is:

```text
.github/workflows/Documentation.yml
```
