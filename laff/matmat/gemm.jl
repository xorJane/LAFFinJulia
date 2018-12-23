"""
    gemm!(alpha::Union{Number, Array{T} where T <: Number}, A::Matrix{T} where T <: Number, B::Matrix{T} where T <: Number, beta::Union{Number, Array{T} where T <: Number}, C::Matrix{T} where T <: Number)

Implementing a blocked, fully optimized matrix-matrix multiply in FlameJulia
would be great, but we haven't made time to do it. We elected
to use Julia's built-in MMM until an ambitious student
takes on the task. 

This routine modifies C with C := alpha * A * B + beta * C
"""

function gemm!(alpha::Union{Number, Array{T} where T <: Number}, A::Matrix{T} where T <: Number, B::Matrix{T} where T <: Number, beta::Union{Number, Array{T} where T <: Number}, C::Matrix{T} where T <: Number)
    # Check parameters
    @assert length(alpha) == 1 "laff.gemm!: alpha must be a number or a 1 element Array"
    @assert length(beta) == 1 "laff.gemm!: beta must be a number or a 1 element Array"
    @assert size(C) == (size(A, 1), size(B, 2)) "laff.gemm!: Dimension mismatch between input and output matrices!"
    @assert size(A, 2) == size(B, 1) "laff.gemm!: Dimension mismatch between input matrices!"
    
	# Update C using the built-in Julia multiply
	C[:,:] = alpha .* A * B + beta .* C
end