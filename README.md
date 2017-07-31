# ILU0.jl

`ILU0.jl` is a Julia implementation of incomplete LU factorization with zero level of fill in. The module is compatible with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

## Requirements

* Julia 0.5 and up

## Instalation

```julia
julia> Pkg.clone("https://github.com/mcovalt/ILU0.jl.git")
```

## Why use ILU0.jl?

You probably shouldn't (thus the unregistered status of the module). Julia's built in factorization methods are much better. Julia uses [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) for sparse matrix factorization which factorizes at about the same speed and results in similarly sized preconditioners which are *much* more robust. In addition, Julia uses heuristics to determine a good factorization scheme for your matrix automatically. In addition, the algorithm used here is for Sparse CSR matrices, therefore an unnecessary transposition must be made. This last point may change in the future.

However, SuiteSparse is a compiled language that does not play nicely with the `Dual` number type found in [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). That's why this was created.

## How to use

```julia
julia> using ILU0
```

* `LU = ilu0(A)`: Create a factorization based on a sparse matrix `A`
* `ilu0!(LU, A)`: Update factorization `LU` in-place based on a sparse matrix `A`. This assumes the original factorization was created with another sparse matrix with the exact same sparsity pattern as `A`. No check is made for this.
* To solve for `x` in `(LU)x=b`, use the same methods as you typically would: `\`, `A_ldiv_B!(x, LU, b)`, `A_ldiv_B!(LU, b)`. See [the docs](https://docs.julialang.org/en/stable/stdlib/linalg/) for further information.

## References

This Julia package follows the algorithm found in [MGMRES](http://people.sc.fsu.edu/~jburkardt/f_src/mgmres/mgmres.html).