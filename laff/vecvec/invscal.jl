"""
    invscal!(alpha, x)

Compute 1/alpha * x, overwriting x

x can be an `Array`, such as a `Vector` or `Matrix`, or a
transposed `Array`. alpha can be a `Number` or an `Array` with
a single numeric element.
"""
function invscal!(alpha::Union{Number, Array{T} where T <: Number}, x::Union{Array{T} where T <: Number, LinearAlgebra.Transpose{T, Array{T}} where T <: Number})
    # Check parameters
    @assert length(alpha) == 1 "laff.invscal!: alpha must be a scalar or an Array with a single element."
    # Get scalar alpha
    alpha = (typeof(alpha) <: Array) ? alpha[1] : alpha
    # Get number elements in x
    n_x = length(x)
    # Update x
    for i in 1:n_x
       x[i] = x[i] / alpha 
    end
end                      