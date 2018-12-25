"""
    trsv!(uplo::String, trans::String, diag::String, A::Matrix{T} where T <: Number, b::Union{LinearAlgebra.Transpose{T, Array{T}} where N where T <: Number, Array{T} where T <:Number} )

Solve A x = b or trans( A ) x = b, overwriting b with x

Parameter uplo indicates whether to use the lower triangular or
upper triangular part of A:
if uplo == "Lower triangular"
   A is lower triangular
elseif upl == "Upper triangular"
   A is upper trianglar

Parameter trans indicates whether to transpose A:
if trans == "No transpose"
   solve A x = b
elseif trans == "Transpose"
   solve trans( A ) x = b

Parameter diag indicates whether A has an (implicit) unit diagonal:
if diag == "Unit diagonal"
   A has an implicit unit diagonal
elseif diag == "Nonunit diagonal"
   Use the entries on the diagonal of A

"""

function trsv!(uplo::String, trans::String, diag::String, A::Matrix{T} where T <: Number, b::Union{LinearAlgebra.Transpose{T,Array{T,N}} where N where T<:Number, Array{T} where T <:Number} )
    # Extract sizes
    m_A, n_A = size(A)
    m_b = length(b)
    
    # Check parameters
    @assert ((uplo == "Lower triangular") || (uplo == "Upper triangular")) "laff.trsv!: illegal value for uplo"
    @assert ((trans == "No transpose") || (trans == "Transpose")) "laff.trsv!: illegal value for trans"
    @assert ((diag == "Nonunit diagonal") || (diag == "Unit diagonal")) "laff.trsv!: illegal value for diag"
    @assert m_b in size(b) "laff.trsv!: b is neither a vector nor a single row or column of a matrix!"
    @assert m_b == n_A "laff.trsv!: size mismatch between b and A ($m_b vs $n_A)"
    
    if uplo == "Lower triangular"
        if trans == "No transpose"
            if diag == "Nonunit diagonal"
                trsv_lnn!( A, b )
            else
                trsv_lnu!( A, b )
            end
        else # trans == "Transpose"
            if diag == "Unit diagonal"
#                 bool = false
#                 (size(b, 2) > size(b, 1)) && (b = Array(transpose(b)); bool = true)
                trsv_ltu!( A, b )
#                 bool && (b = Array(transpose(b)))
            else
                println("laff.trsv!: trans == Transpose not yet implemented for Lower triangular, with Nonunit diagonal")
            end
        end
        
    else # uplo == "Upper triangular"
        if trans == "No transpose"
            if diag == "Nonunit diagonal"
                trsv_unn!( A, b )
            else
                trsv_unu!( A, b )
            end
        else # trans == "Transpose"
            if diag == "Nonunit diagonal"
                trsv_utn!( A, b )
#                 bool = false
#                 (size(b, 2) > size(b, 1)) && (b = Array(transpose(b)); bool = true)
#                 trsv_utn!( A, b )
#                 bool && (b = Array(transpose(b)))
            else
                println("laff.trsv!: trans == Transpose not yet implemented for Upper Triangular with Unit diagonal")
            end
        end        
    end
end
