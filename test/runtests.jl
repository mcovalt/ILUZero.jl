#!/usr/bin/env julia

#Start Test Script
using ILU0
using BiCGStab
using Base.Test

type BiCGStabFunc{T<:Real} <: AbstractBiCGStabFunc{T}
    A::SparseMatrixCSC{T,Int64}
end
function (t::BiCGStabFunc){T<:Real}(out::Vector{T}, x::Vector{T})
    A_mul_B!(out, t.A, x)
end
type BiCGStabNoPrecon{T<:Real} <: AbstractBiCGStabPrecon{T} end
function (t::BiCGStabNoPrecon){T<:Real}(out::Vector{T}, b::Vector{T})
    out .= b
end
type BiCGStabPrecon{T<:Real} <: AbstractBiCGStabPrecon{T}
    LU::ILU0Precon{T,Int64}
end
function (t::BiCGStabPrecon){T<:Real}(out::Vector{T}, b::Vector{T})
    A_ldiv_B!(out, t.LU, b)
end

function test_solve()
    A = sprandn(100,100,.1) + 10.0*speye(100)
    BCGA = BiCGStabFunc(A)
    BCGNM = BiCGStabNoPrecon{Float64}()
    BCGM = BiCGStabPrecon(ilu0(A))
    BCGD = BiCGStabData(100, Float64)
    x = zeros(100)
    b = rand(100)
    true_x = A\b
    outint, nocon_niter = BiCGStab!(BCGA, BCGNM, BCGD, x, b)
    println("No preconditioning: ", nocon_niter, " iterations")
    x .= 0.0
    outint, con_niter = BiCGStab!(BCGA, BCGM, BCGD, x, b)
    println("Preconditioned: ", con_niter, " iterations")
    if nocon_niter > con_niter
        return true
    else
        return false
    end
end

@test test_solve()