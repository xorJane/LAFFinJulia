"""
    zerov!(x::Union{Matrix{T} where T <:Number, Vector{T} where T <:Number, LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, LinearAlgebra.Transpose{T, Matrix{T}} where T <: Number})

Set all components of x to zero.

x contains elements of type `Number` and can be a `Vector`, a transposed `Vector` or a (transposed) `Matrix` with at least one unary dimension.
"""
function zerov!(x::Union{Matrix{T} where T <:Number, Vector{T} where T <:Number, LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, LinearAlgebra.Transpose{T, Matrix{T}} where T <: Number})
    xsize = size(x)
    length(xsize) > 1 ? ( @assert (1 in xsize) "Cannot pass Matrix with both multiple columns and rows to laff.zerov!") : 0
    n = length(x)
    # Assign elements of x to zero of the same type
    # of the elements of the original x
    correctzero = (n > 0) ? zero(x[1]) : 0
    for i in 1:n
        x[i] = correctzero
    end
end
