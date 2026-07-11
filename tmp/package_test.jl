using Pkg

root = "/home/jgkong/project/2026/0710/register-grassmann-package"
Pkg.activate(root)
Pkg.instantiate()
Pkg.test()

write(joinpath(root, "success.sentinel"), "grassmann_symbolics_pkg_test_passed\n")
