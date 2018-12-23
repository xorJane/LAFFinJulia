"""
    trsm!(uplo::String, trans::String, diag::String, A::Matrix{T} where T <: Number, B::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number})

Solve A X = B or trans( A X ) = trans( B ), overwriting B with X

Parameter uplo indicates whether to use the lower triangular or
upper triangular part of A:
if uplo == "Lower triangular"
   A is lower triangular
elseif uplo == "Upper triangular"
   A is upper trianglar

Parameter trans indicates whether to transpose A:
if trans == "No transpose"
   solve A X = B
elseif trans == "Transpose"
   solve trans( A X ) = trans( B )

Parameter diag indicates whether A has an (implicit) unit diagonal:
if diag == "Unit diagonal"
   A has an implicit unit diagonal
elseif diag == "Nonunit diagonal"
   Use the entries on the diagonal of A

"""
function trsm!(uplo::String, trans::String, diag::String, A::Matrix{T} where T <: Number, B::Union{LinearAlgebra.Transpose{T, Matrix{T}} where T <: Number, Matrix{T} where T <:Number})
    # Extract sizes
    m_A, n_A = size(A)
    m_B, n_B = size(B)
    
    # Check parameters
    @assert ((uplo == "Lower triangular") || (uplo == "Upper triangular")) "laff.trsm!: illegal value for uplo"
    @assert ((trans == "No transpose") || (trans == "Transpose")) "laff.trsm!: illegal value for trans"
    @assert ((diag == "Nonunit diagonal") || (diag == "Unit diagonal")) "laff.trsm!: illegal value for diag"
    @assert m_B == n_A "laff.trsm!: size mismatch between B and A"
    
    if uplo == "Lower triangular"
        if trans == "No transpose"
            if diag == "Unit diagonal"
                trsm_lnu!( A, B )
            else
                println( "laff.trsm!: diag == Nonunit diagonal not yet implemented for Lower triangular" )
            end
        else # trans == "Transpose"
            if diag == "Unit diagonal"
                trsm_ltu!( A, B )
            else
                println( "laff.trsm!: trans == Transpose not yet implemented for Lower triangular, nonunit diagonal" )
            end
        end
        
    else # uplo == "Upper triangular"
        if trans == "No transpose"
            if diag == "Unit diagonal"
                println( "laff.trsm!: trans == No transpose not yet implemented for Upper triangular, unit diagonal" )
            else
                trsm_unn!( A, B )
            end
        else # trans == "Transpose"
            if diag == "Nonunit diagonal"
                trsm_utn!( A, B )
            else
                println( "laff.trsm: diag == Unit diagonal not yet implemented for Upper triangular" )
            end
        end        
    end
end
