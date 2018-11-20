"""
    gemv!(trans, alpha, A, x, beta, y)

Compute y := alpha * A * x + beta * y

x and y can be `Array`s or transposed `Array`s. For example, x and y can be column vectors (`Vector`s) or row vectors (transposed `Vector`s).   If necessary, a
transposition happens with those Arrays. A is a `Matrix` or transposed `Matrix. `alpha` and `beta` are numbers or matrices with single numeric elements.
"""

function gemv!(trans::String, alpha::Union{Number, Matrix{T} where T <: Number}, A::Union{Matrix{T} where T <:Number, LinearAlgebra.Transpose{T, Matrix{T}} where T <: Number}, x::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number}, beta::Union{Number, Matrix{T} where T <: Number}, y::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number})
    
    # Check parameters
    @assert (trans == "No transpose" || trans == "Transpose") "laff.gemv!: illegal value for trans"
    @assert (length(alpha) == 1) "laff.gemv!: illegal value for alpha"
    @assert (length(beta) == 1) "laff.gemv!: illegal value for beta"
    
    # Grab scalar versions of alpha and beta
    α = (typeof(alpha) <: Matrix) ? alpha[1, 1] : alpha
    β = (typeof(beta) <: Matrix) ? beta[1, 1] : beta
    
    # Transpose x if it is either a transposed object or a matrix with a single row.
    # Same for y.
    x = (typeof(x) <: LinearAlgebra.Transpose || (typeof(x) <: Matrix && size(x, 1) == 1)) ? transpose(x) : x
    y = (typeof(y) <: LinearAlgebra.Transpose || (typeof(y) <: Matrix && size(y, 1) == 1)) ? transpose(y) : y
    
    # Extract sizes         
    n_x = length(x) 
    n_y = length(y)
    m_A, n_A = size(A)

    # Update y
    if "No transpose" == trans
        @assert m_A == n_y "laff.gemv!: size mismatch between y and A"
        @assert n_A == n_x "laff.gemv!: size mismatch between x and A"
        y .= α * A * x + β * y

    else
        @assert m_A == n_x "laff.gemv!: size mismatch between x and A"
        @assert n_A == n_y "laff.gemv!: size mismatch between y and A"
        y .= α * transpose(A) * x + β * y
    end              
end
