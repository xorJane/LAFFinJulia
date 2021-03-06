"""
    dot( x::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number}, y::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number}, alpha::Number = 0.0)

Compute alpha = x^T + y, storing the result in alpha and returning alpha.

x and y can be `Array`s or transposed `Array`s. For example, x and y can be column vectors (`Vector`s) or row vectors (transposed `Vector`s). If necessary, an implicit transposition happens.
"""
function dot(x::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number}, y::Union{LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, Vector{T} where T <:Number}, alpha::Number = 0.0)
    m = length(x)
    @assert m == length(y) "Dimension mismatch between x and y."
    for i in 1:m;  alpha += x[i] * y[i]; end
    return alpha
end
