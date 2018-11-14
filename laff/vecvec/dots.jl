"""
    dot!(x, y, alpha)

Compute alpha = x' * y + alpha, storing the result in alpha and updating alpha.

x and y can be `Array`s or transposed `Array`s. For example, x and y can be column vectors (`Vector`s) or row vectors (transposed `Vector`s). If necessary, an implicit transposition happens. alpha is an `Array` containing a single number.
"""
function dot!(x::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number}, y::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number}, alpha::Array{T} where T <: Number)

    @assert length(alpha) == 1 "laff.dot!: alpha must contain one element."

    m = length(x)
    @assert m == length(y) "Dimension mismatch between x and y."

    for i in 1:m
        alpha[1] += x[i] * y[i]
    end
end
