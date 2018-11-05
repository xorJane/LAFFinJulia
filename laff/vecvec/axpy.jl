using LinearAlgebra
"""
    axpy!(alpha::Number, x::Union{LinearAlgebra.Transpose{Vector}, Vector}, y::Union{LinearAlgebra.Transpose{Vector}, Vector})

Compute y = alpha * x + y, overwriting y.

x and y can be row and/or column vectors. If necessary, an implicit transposition happens.
"""
function axpy!(alpha::Number, x::Union{LinearAlgebra.Transpose{Vector}, Vector}, y::Union{LinearAlgebra.Transpose{Vector}, Vector})
    m = length(x)
    @assert m == length(y) "Dimension mismatch between x and y."
    for i in 1:m; y[i] += alpha * x[i]; end
end


