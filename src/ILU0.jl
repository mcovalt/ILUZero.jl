module ILU0
    # Overloaded functions
    import Base: \, A_ldiv_B!

    # ILU0 type definition
    immutable ILU0Precon{T<:Real, N<:Integer} <: Factorization{T}
        n::Int64
        mat::SparseMatrixCSC{T,N}
        rowdiag::Vector{N}
        iw::Vector{N}
    end

    # Allocates ILU0Precon type
    function allocILU0Precon{T<:Real, N<:Integer}(A::SparseMatrixCSC{T,N})
        mat = A'
        rowptr = mat.colptr
        colval = mat.rowval
        nzval = mat.nzval
        n = size(A)[1]
        rowdiag = zeros(N, n)
        iw = zeros(N, n)
        diagFound = false
        for i = 1:n
            for k = rowptr[i]:rowptr[i+1] - 1
                if colval[k] == i
                    rowdiag[i] = k
                    diagFound = true
                    break
                end
            end
            if diagFound == false
                error("Diagonal element missing in sparse matrix.")
            else
                diagFound = false
            end
        end
        return ILU0Precon(n, mat, rowdiag, iw)
    end

    # Updates ILU0Precon type in-place based on matrix A
    function ilu0!{T<:Real, N<:Integer}(LU::ILU0Precon{T,N}, A::SparseMatrixCSC{T,N})
        transpose!(LU.mat, A)
        rowptr = LU.mat.colptr
        colval = LU.mat.rowval
        nzval = LU.mat.nzval
        j = 0
        i = 0
        jrow = 0
        for i = 1:LU.n
            LU.iw .= -1
            for k = rowptr[i]:rowptr[i+1] - 1
                LU.iw[colval[k]] = k
            end
            for j = rowptr[i]:rowptr[i+1] - 1
                jrow = colval[j]
                if i <= jrow
                    break
                end
                tl = nzval[j]*nzval[LU.rowdiag[jrow]]
                nzval[j] = tl
                for jj = LU.rowdiag[jrow] + 1:rowptr[jrow+1] - 1
                    jw = LU.iw[colval[jj]]
                    if jw != -1
                        nzval[jw] -= tl*nzval[jj]
                    end
                end
            end
            if jrow != i
                error("Incompatible sparse matrix (jrow != i")
            elseif nzval[j] == 0.0
                error("Zero pivot.")
            end
            nzval[j] ^= -1.0
        end
        nzval[LU.rowdiag] .^= -1.0
        return
    end

    # Constructs ILU0Precon type based on matrix A
    function ilu0{T<:Real, N<:Integer}(A::SparseMatrixCSC{T,N})
        size(A)[1] == size(A)[2] || error("Sparse matrix must be square.")
        LU = allocILU0Precon(A)
        ilu0!(LU, A)
        return LU
    end

    # Solves LU\b overwriting b
    function A_ldiv_B!{T<:Real, N<:Integer}(LU::ILU0Precon{T,N}, b::Vector{T})
        length(b) == LU.n || throw(DimensionMismatch())
        rowptr = LU.mat.colptr
        colval = LU.mat.rowval
        nzval = LU.mat.nzval
        for i = 2:LU.n
            for j = rowptr[i]:LU.rowdiag[i] - 1
                b[i] -= nzval[j]*b[colval[j]]
            end
        end
        for i = LU.n:-1:1
            for j = LU.rowdiag[i] + 1:rowptr[i+1] - 1
                b[i] -= nzval[j]*b[colval[j]]
            end
            b[i] /= nzval[LU.rowdiag[i]]
        end
        return
    end

    # Solves LU\b overwriting x
    function A_ldiv_B!{T<:Real, N<:Integer}(x::Vector{T}, LU::ILU0Precon{T,N}, b::Vector{T})
        (length(x) == LU.n && length(b) == LU.n) || throw(DimensionMismatch())
        x .= b
        Base.A_ldiv_B!(LU, x)
        return
    end

    # Returns LU\b
    function \{T<:Real, N<:Integer}(LU::ILU0Precon{T,N}, b::Vector{T})
        length(b) == LU.n || throw(DimensionMismatch())
        x = zeros(T, length(b))
        A_ldiv_B!(x, LU, b)
        return x
    end

    export ILU0Precon, \, A_ldiv_B!, ilu0, ilu0!
end
