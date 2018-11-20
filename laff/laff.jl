module laff
    using LinearAlgebra
    include("vecvec/copy.jl")
    include("vecvec/scal.jl")
    include("vecvec/axpy.jl")
    include("vecvec/dot.jl")
    include("vecvec/norm2.jl")
    include("util/zerov.jl")
    include("util/onev.jl")
    include("vecvec/dots.jl")
    include("matvec/gemv.jl")
    include("matvec/ger.jl")

end
