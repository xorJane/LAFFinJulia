"""
    copy!(x, y)

Compute y = x, overwriting y.

x and y can be arrays with the same number of elements. For example,
x and y may both be columns vectors (`Vector`s) and/or row vectors
(transposed `Vector`s) of the same length. 
"""
function copy!( x, y )
    n = length(x)
    @assert n == length(y) "x and y have different numbers of elements!"
    
    for i in 1:n
        y[ i ] = x[ i ]
    end
end
