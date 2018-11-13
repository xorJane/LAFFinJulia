"""
    onev!(x::Union{Matrix{T} where T <:Number, Vector{T} where T <:Number, LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, LinearAlgebra.Transpose{T, Matrix{T}} where T <: Number})

Set all components of x to one.

x contains elements of type `Number` and can be a `Vector`, a transposed `Vector` or a (transposed) `Matrix` with at least one unary dimension.
"""
function onev!(x::Union{Matrix{T} where T <:Number, Vector{T} where T <:Number, LinearAlgebra.Transpose{T, Vector{T}} where T <: Number, LinearAlgebra.Transpose{T, Matrix{T}} where T <: Number})
    xsize = size(x)
    length(xsize) > 1 ? ( @assert (1 in xsize) "Cannot pass Matrix with both multiple columns and rows to laff.onev!") : 0
    n = length(x)
    # Assign elements of x to one of the same type
    # of the elements of the original x
    correctone = (n > 0) ? one(x[1]) : 1
    for i in 1:n
        x[i] = correctone
    end
end
