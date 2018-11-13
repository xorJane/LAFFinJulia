"""
    copy!( x::Union{Array{T} where T <:Number, LinearAlgebra.Transpose{T, Array{T}} where T <: Number}, y::Union{Array{T} where T <:Number, LinearAlgebra.Transpose{T, Array{T}} where T <: Number} )

Compute y = x, overwriting y.

x and y can be `Array`s or transposed `Array`s with the same number of elements. For example,
x and y may both be columns vectors (`Vector`s) and/or row vectors
(transposed `Vector`s) of the same length.
"""
function copy!( x::Union{Array{T} where T <:Number, LinearAlgebra.Transpose{T, Array{T}} where T <: Number}, y::Union{Array{T} where T <:Number, LinearAlgebra.Transpose{T, Array{T}} where T <: Number} )
    n = length(x)
    @assert n == length(y) "x and y have different numbers of elements!"

    for i in 1:n
        y[ i ] = x[ i ]
    end
end
