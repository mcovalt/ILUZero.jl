#!/usr/bin/env julia

#Start Test Script
using ILU0, IterativeSolvers, LinearAlgebra, SparseArrays, Test

function test_solve()
    n = 100
    tA = sprandn(n,n,.1) + 10.0*I
    A = tA'*tA
    ilu_prec = ilu0(A)
    b = rand(n)
    x, ch = cg(A, b, log=true)
    nocon_niter = ch.iters
    println("No preconditioning: ", nocon_niter, " iterations")
    x, ch = cg(A, b, Pl=ilu_prec, log=true)
    con_niter = ch.iters
    println("Preconditioned: ", con_niter, " iterations")
    if nocon_niter > con_niter
        return true
    else
        return false
    end
    return true
end

@test test_solve()
