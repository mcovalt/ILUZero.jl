#!/usr/bin/env julia

#Start Test Script
using ILU0
using IterativeSolvers
using Base.Test

function test_solve()
    tA = sprandn(100,100,.1) + 10.0*speye(100)
    A = tA'*tA
    LU = ilu0(A)
    b = rand(100)
    x, ch = cg(A, b, log=true)
    nocon_niter = ch.iters
    println("No preconditioning: ", nocon_niter, " iterations")
    x, ch = cg(A, b, Pl=LU, log=true)
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