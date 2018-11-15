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

"""
    scal!(alpha, x)

Compute alpha * x, overwriting x.

alpha must be a `Matrix` with one numeric element. x can be an `Array` or the transpose of an Array.
"""
function scal!(alpha::Matrix{T} where T <:Number, x::Union{Array{T} where T <:Number, LinearAlgebra.Transpose{T, Array{T}} where T <: Number})
    @assert length(alpha) == 1 "laff.scal!: alpha should be a Number or a Matrix with one numeric element!"
    m = length(x)
    for i in 1:m
        x[i] = alpha[1] * x[i]
    end
end
