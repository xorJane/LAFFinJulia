using LinearAlgebra
"""
    dot( x::Union{LinearAlgebra.Transpose{Vector}, Vector}, y::Union{LinearAlgebra.Transpose{Vector}, Vector}, alpha::Number = 0.0)

Compute alpha = x^T + y, storing the result in alpha and returning alpha.

x and y can be row and/or column vectors. If necessary, an implicit transposition happens.
"""
function dot(x::Union{LinearAlgebra.Transpose{Vector}, Vector}, y::Union{LinearAlgebra.Transpose{Vector}, Vector}, alpha::Number = 0.0)
    m = length(x)
    @assert m == length(y) "Dimension mismatch between x and y."
    for i in 1:m;  alpha += x[i] * y[i]; end
    return alpha
end


