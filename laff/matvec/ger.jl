"""
    ger!(alpha, y, x, A)

Compute A := alpha * y * x + A

x and y can be row and/or column vectors.  If necessary, a
transposition happens.
"""
function ger!(alpha::Union{Number, Matrix{T} where T <: Number}, y::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number}, x::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number}, A::Union{Matrix{T} where T <:Number, LinearAlgebra.Transpose{T, Matrix{T}} where T <: Number})
    # Check parameters
    @assert (length(alpha) == 1) "laff.gemv!: illegal value for alpha"
    
    # Grab scalar version of alpha
    α = (typeof(alpha) <: Matrix) ? alpha[1, 1] : alpha
    
    # Make sure x is a row of a matrix or a transposed Vector
    # and that y is a column of a matrix or a Vector
    # Transpose either or both as necessary.
    x = (typeof(x) <: LinearAlgebra.Transpose || (typeof(x) <: Matrix && size(x, 1) == 1)) ? x : transpose(x)
    y = (typeof(y) <: LinearAlgebra.Transpose || (typeof(y) <: Matrix && size(y, 1) == 1)) ? transpose(y) : y
    
    # Extract sizes
    n_x = length(x) 
    n_y = length(y)
    m_A, n_A = size(A)
   
    # Update A
    A .= α * y * x + A
end