using LinearAlgebra

"""
    scal!(alpha, x)

Compute alpha * x, overwriting x.

alpha must be a numnber and x can be row or column vectors.
"""
function scal!(alpha::Number, x::Union{LinearAlgebra.Transpose, Array})
    m = length(x)
    for i in 1:m
        x[i] = alpha * x[i]
    end
end
