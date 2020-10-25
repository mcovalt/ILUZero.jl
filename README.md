# ILUZero.jl

`ILUZero.jl` is a Julia implementation of incomplete LU factorization with zero level of fill-in. It allows for non-allocating updates of the factorization. The module is compatible with [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

## Requirements

* Julia 1.0 and up

## Installation

```julia
julia> ]
pkg> add ILUZero
```

## Why use ILUZero.jl?

You probably shouldn't. Julia's built in factorization methods are much better. Julia uses [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) for sparse matrix factorization which factorizes at about nearly the same speed and results in similarly sized preconditioners which are *much* more robust. In addition, Julia uses heuristics to determine a good factorization scheme for your matrix automatically.

Due to the zero-fill of this package, however, factorization should be a bit faster and preconditioners can be preallocated if updated by a matrix of identical sparsity.

## How to use

```julia
julia> using ILUZero
```

* `LU = ilu0(A)`: Create a factorization based on a sparse matrix `A`
* `ilu0!(LU, A)`: Update factorization `LU` in-place based on a sparse matrix `A`. This assumes the original factorization was created with another sparse matrix with the exact same sparsity pattern as `A`. No check is made for this.
* To solve for `x` in `(LU)x=b`, use the same methods as you typically would: `\` or `ldiv!(x, LU, b)`. See [the docs](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) for further information.
* There's also:
  - Forward substitution: `forward_substitution!(y, LU, b)` solves `L\b` and stores the solution in y.
  - Backward substitution: `backward_substitution!(x, LU, y)` solves `U\y` and stores the solution in x.
  - Nonzero count: `nnz(LU)` will return the number of nonzero entries in `LU`.

## Performance

```julia
julia> using ILUZero
julia> using BenchmarkTools, LinearAlgebra, SparseArrays
julia> A = sprand(1000, 1000, 5 / 1000) + 10I
julia> fact = @btime ilu0(A)
       107.600 μs (16 allocations: 160.81 KiB)
julia> updated_fact = @btime ilu0!($fact, $A)
       71.500 μs (0 allocations: 0 bytes)
```
