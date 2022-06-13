# ILUZero.jl

`ILUZero.jl` is a Julia implementation of incomplete LU factorization with zero level of fill-in. It allows for non-allocating updates of the factorization.

## Requirements

* Julia 1.0 and up

## Installation

```julia
julia> ]
pkg> add ILUZero
```

## Why use ILUZero.jl?

This package was created to exploit some specific properties of a problem I had:
- the non-zero elements of the matrix were frequently updated
- an approximation of the matrix inverse was needed each time the non-zero elements were changed
- the sparsity structure of the matrix was never altered
- the non-zero elements of the matrix weren't necessarily floating point numbers

The ILU(0) preconditioner was a great solution for this problem. The method is simple enough to write a performant version purely in Julia. This allowed using Julia's flexible type system. Also, ILU(0) has constant memory requirements so the updating matrix could re-use previously allocated structures<sup>*</sup>. Win win (win win)!

_<sup>*</sup>Hint: take a look at [ConjugateGradients.jl](https://github.com/mcovalt/ConjugateGradients.jl) if this type of preconditioner would be beneficial to your project. Combined, these packages can help reduce allocations in those hot paths._



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
