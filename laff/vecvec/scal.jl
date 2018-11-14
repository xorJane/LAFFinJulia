"""
    scal!(alpha, x)

Compute alpha * x, overwriting x.

alpha must be a `Number`. x can be an `Array` or the transpose of an Array.
"""
function scal!(alpha::Number, x::Union{Array{T} where T <:Number, LinearAlgebra.Transpose{T, Array{T}} where T <: Number})
    m = length(x)
    for i in 1:m
        x[i] = alpha * x[i]
    end
end
