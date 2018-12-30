"""
    norm2( x::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number} )

Compute the 2-norm of an Array, returning alpha. x can be an `Array` or a transposed `Array`.
"""
function norm2( x::Union{LinearAlgebra.Transpose{T, Array{T}} where T <: Number, Array{T} where T <:Number} )
    # Ensure that we don't modify x in any way by copying it to a new vector, y
    m = length(x)
    y = fill(0.0, m)
    laff.copy!(x, y)

    # Initialize vairables that we will use to appropriate values
    alpha = 0.0
    maxval = y[1]

    # Find a value to scale by to avoid under/overflow
    for i in 1:m
        if abs(y[i]) > maxval
            maxval = abs(y[i])
        end
    end

    # If y is the zero vector, return 0.0
    if abs(maxval) < 10^-7
        return 0
    end

    # Scale all of the values by 1/maxval to prevent under/overflow
    laff.scal!(1.0/maxval, y)

    alpha = maxval * sqrt(laff.dot(y, y))

    return alpha
end
