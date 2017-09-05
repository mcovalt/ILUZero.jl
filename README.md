# ILU0.jl

`ILU0.jl` is a Julia implementation of incomplete LU factorization with zero level of fill-in. It allows for non-allocating updates of the factorization. The module is compatible with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

## Requirements

* Julia 0.6 and up

## Instalation

```julia
julia> Pkg.clone("https://github.com/mcovalt/ILU0.jl.git")
```

## Why use ILU0.jl?

You probably shouldn't (thus the unregistered status of the module). Julia's built in factorization methods are much better. Julia uses [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) for sparse matrix factorization which factorizes at about nearly the same speed and results in similarly sized preconditioners which are *much* more robust. In addition, Julia uses heuristics to determine a good factorization scheme for your matrix automatically.

However, SuiteSparse is a compiled language that does not play nicely with the `Dual` number type found in [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). That's why this was created.

The package [ILU.jl](https://github.com/haampie/ILU.jl) is also written in pure Julia and ought to be compatible with ForwardDiff.jl as well. It allows a drop tolerance to be set so the fill-in of the preconditioner can be above zero which ought to result in a more robust preconditioner. Due to the zero-fill of this package, however, factorization should be a bit faster and preconditioners can be preallocated if updated by a matrix of identical sparsity.

## How to use

```julia
julia> using ILU0
```

* `LU = ilu0(A)`: Create a factorization based on a sparse matrix `A`
* `ilu0!(LU, A)`: Update factorization `LU` in-place based on a sparse matrix `A`. This assumes the original factorization was created with another sparse matrix with the exact same sparsity pattern as `A`. No check is made for this.
* To solve for `x` in `(LU)x=b`, use the same methods as you typically would: `\` or `A_ldiv_B!(x, LU, b)`. See [the docs](https://docs.julialang.org/en/stable/stdlib/linalg/) for further information.

## Performance

```julia
julia> using ILU0
julia> using BenchmarkTools
julia> A = sprand(1000, 1000, 5 / 1000) + 10I
julia> fact = @btime ilu0(A)
       158.973 μs (16 allocations: 164.94 KiB)
julia> updated_fact = @btime ilu0!($fact, $A)
       105.015 μs (0 allocations: 0 bytes)
```
