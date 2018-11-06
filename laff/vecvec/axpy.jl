"""
    axpy!(alpha::Number, x::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number}, y::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number})

Compute y = alpha * x + y, overwriting y.

x and y can be row and/or column vectors. If necessary, an implicit transposition happens.
"""
function axpy!(alpha::Number, x::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number}, y::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number})
    m = length(x)
    @assert m == length(y) "Dimension mismatch between x and y."
    for i in 1:m; y[i] += alpha * x[i]; end
end
