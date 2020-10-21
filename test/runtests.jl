# Start Test Script
using ILU0, LinearAlgebra, SparseArrays, Test

function cg(A, b; M=I)
    n = size(A, 1)
    x = zeros(n)
    r = copy(b)
    z = zeros(n)
    ldiv!(z, M, r)
    p = copy(z)
    γ = dot(r, z)
    k = 0
    tired = false
    solved = false
    while !(solved || tired)
        Ap = A * p
        α = γ / dot(p, Ap)
        x = x + α * p
        r = r - α * Ap
        ldiv!(z, M, r)
        γ_new = dot(r, z)
        β = γ_new / γ
        γ = γ_new
        p = z + β * p
        k = k + 1
        tired = k > n
        solved = norm(r) ≤ 1e-8 * norm(b)
    end
    return x, k
end

function test_solve()
    n = 100
    tA = sprandn(n, n, .1) + 10.0 * I
    A = tA' * tA
    ilu_prec = ilu0(A)
    b = rand(n)
    x1, k1 = cg(A, b)
    println("No preconditioning: ", k1, " iterations")
    x2, k2 = cg(A, b, M=ilu_prec)
    println("Preconditioned: ", k2, " iterations")
    return k1 > k2
end

function test_substitutions()
    n = 100
    tA = sprandn(n, n, .1) + 10.0 * I
    A = tA' * tA
    ilu_prec = ilu0(A)
    b = rand(n)
    x = ilu_prec \ b
    x1 = zeros(n)
    x2 = zeros(n)
    forward_substitution(x1, ilu_prec, b)
    backward_substitution(x2, ilu_prec, x1)
    return (x2 == x)
end

@test test_solve()
@test test_substitutions()
