module flame

using LinearAlgebra

"""
    merge_2x2!(TL::Matrix, TR::Matrix, BL::Matrix, BR::Matrix, A::Matrix)

Modify `A` to contain the values of submatrices/quadrants TL, TR,
BL, and BR.
"""
function merge_2x2!(TL::Matrix, TR::Matrix, BL::Matrix, BR::Matrix, A::Matrix)
    if size(TL)[1] > 0 && size(TL)[2] > 0
        for i in 1:size(TL)[1], j in 1:size(TL)[2]
           A[i, j] = TL[i, j]
        end
    end

    if size(TR)[1] > 0 && size(TR)[2] > 0
        for i in 1:size(TR)[1], j in 1:size(TR)[2]
           A[i, j + size(TL)[2]] = TR[i, j]
        end
    end

    if size(BL)[1] > 0 && size(BL)[2] > 0
        for i in 1:size(BL)[1], j in 1:size(BL)[2]
           A[i + size(TL)[1], j] = BL[i, j]
        end
    end

    if size(BR)[1] > 0 && size(BR)[2] > 0
        for i in 1:size(BR)[1], j in 1:size(BR)[2]
           A[i + size(TL)[1], j + size(TL)[2]] = BR[i, j]
        end
    end
end

"""
    merge_2x1!(T::Vector, B::Vector, x::Vector)

Combine input vectors `T` and `B` to create input/output vector `x`.
"""
function merge_2x1!(T::Vector, B::Vector, x::Vector)
    xtemp = vcat(T, B)
    for i in 1:length(xtemp)
        x[i] = xtemp[i]
    end
end

"""
    merge_2x1!(T::Matrix, B::Matrix, A::Union{Matrix, LinearAlgebra.Transpose{T, Matrix{T}} where T})

Combine input matrices `T` and `B` vertically to create input/output matrix `A`.
"""
function merge_2x1!(T::Matrix, B::Matrix, A::Union{Matrix, LinearAlgebra.Transpose{T, Matrix{T}} where T})
    Atemp = vcat(T, B)
    m, n = size(Atemp)
    for i in 1:m, j in 1:n
        A[i, j] = Atemp[i, j]
    end
end

"""
    merge_1x2!(L::Matrix, R::Matrix, A::Matrix)

Update matrix `A` to contain the elements of the submatrices
`L` on the left and `R` on the right.
"""
function merge_1x2!(L::Matrix, R::Matrix, A::Matrix)
    if size(L)[1] > 0 && size(L)[2] > 0
        for i in 1:size(L)[1], j in 1:size(L)[2]
            A[i, j] = L[i, j]
        end
    end

    if size(R)[1] > 0 && size(R)[2] > 0
        for i in 1:size(R)[1], j in 1:size(R)[2]
           A[i, j + size(L)[2]] = R[i, j]
        end
    end
end

"""
    cont_with_1x3_to_1x2(A1::Matrix, A2::Matrix, A3::Matrix, side = "LEFT")

Repartition three matrices (vertical slabs of an original matrix) into two by combining the middle portion, A2, of the original matrix with the `side` input submatrix.
"""
function cont_with_1x3_to_1x2(A1::Matrix, A2::Matrix, A3::Matrix, side = "LEFT")
    (side == "LEFT") ? (hcat(A1, A2), A3) : (A1, hcat(A2, A3))
end

"""
    cont_with_3x1_to_2x1(A1::Vector, A2::Vector, A3::Vector, side = "TOP")

Repartition three vectors into two by combining the middle portion, A2, of the original vector with the `side` input subvector.

"""
function cont_with_3x1_to_2x1(A1::Vector, A2::Vector, A3::Vector, side = "TOP")
    (side == "TOP") ? (vcat(A1, A2), A3) : (A1, vcat(A2, A3))
end

"""
    cont_with_3x1_to_2x1(A1::Matrix, A2::Matrix, A3::Matrix, side = "TOP")

Repartition three matrices (horizontal slabs of an original matrix) into two by combining the middle portion, A2, of the original matrix with the `side` input submatrix.
"""
function cont_with_3x1_to_2x1(A1::Matrix, A2::Matrix, A3::Matrix, side = "TOP")
    (side == "TOP") ? (vcat(A1, A2), A3) : (A1, vcat(A2, A3))
end

"""
    cont_with_3x3_to_2x2(A11::Matrix, A12::Matrix, A13::Matrix,
                         A21::Matrix, A22::Matrix, A23::Matrix,
                         A31::Matrix, A32::Matrix, A33::Matrix, quad = "TL")

Concatenate matrices together to repartition the original matrix into 4 quadrants
rather than a 3x3 grid. The middle submatrix `A22` is included in the `quad` output quadrant.
"""
function cont_with_3x3_to_2x2(A11::Matrix, A12::Matrix, A13::Matrix,
                              A21::Matrix, A22::Matrix, A23::Matrix,
                              A31::Matrix, A32::Matrix, A33::Matrix, quad = "TL")
    if quad == "TL"
        TL = hcat(vcat(A11, A21), vcat(A12, A22))
        TR = vcat(A13, A23)
        BL = hcat(A31, A32)
        BR = A33
    elseif quad == "TR"
        TL = vcat(A11, A21)
        TR = hcat(vcat(A12, A22), vcat(A13, A23))
        BL = A31
        BR = hcat(A32, A33)
    elseif quad == "BL"
        TL = hcat(A11, A12)
        TR = A13
        BL = hcat(vcat(A21, A31), vcat(A22, A32))
        BR = vcat(A23, A33)
    elseif quad == "BR"
        TL = A11
        TR = hcat(A12, A13)
        BL = vcat(A21, A31)
        BR = hcat(vcat(A22, A32), vcat(A23, A33))
    end

    return TL, TR, BL, BR
end

"""
    part_1x2(A::Matrix, vpart = 0, side = "LEFT")

Partition a matrix into left and righthand portions, with
`vpart` columns in the `side`hand side.
"""
function part_1x2(A::Matrix, vpart = 0, side = "LEFT")
    if vpart < 0
        throw(DimensionMismatch("size < 0"))
    elseif vpart > size(A)[2]
        throw(DimensionMismatch("size > col dimension"))
    elseif !(side in ("LEFT", "RIGHT"))
        throw(ArgumentError("""side must be "LEFT" or "RIGHT" """))
    end

    vpart = (side == "LEFT") ? vpart : size(A)[2] - vpart

    AL, AR = A[:, 1:vpart], A[:, vpart + 1:end]
end

"""
    part_2x1(x::Vector, hpart=0, side="TOP")

Partition a vector into top and bottom portions, with
`hpart` elements in the `side` portion.
"""
function part_2x1(x::Vector, hpart=0, side="TOP")
    if hpart < 0
        throw(DimensionMismatch("size < 0"))
    elseif hpart > size(x)[1]
        throw(DimensionMismatch("size > row dimension"))
    elseif !(side in ("TOP", "BOTTOM"))
        throw(ArgumentError("""side must be "TOP" or "BOTTOM" """))
    end
    hpart = (side == "TOP") ? hpart : size(x)[1] - hpart
    xT, xB = x[1:hpart], x[hpart + 1:end]
end

"""
    part_2x1(A::Matrix, hpart = 0, side = "TOP")

Partition a matrix into top and bottom portions, with
`hpart` rows in the `side` portion.
"""
    function part_2x1(A::Union{Matrix, LinearAlgebra.Transpose{T, Matrix{T}} where T}, hpart=0, side="TOP")
    if hpart < 0
        throw(DimensionMismatch("size < 0"))
    elseif hpart > size(A)[1]
        throw(DimensionMismatch("size > row dimension"))
    elseif !(side in ("TOP", "BOTTOM"))
        throw(ArgumentError("""side must be "TOP" or "BOTTOM" """))
    end
    hpart = (side == "TOP") ? hpart : size(A)[1] - hpart
    AT, AB = A[1:hpart, :], A[hpart + 1:end, :]
end

"""
    part_2x2(A::Matrix, m::Int, n::Int, quad::String)

Break `A` into four quadrants. `quad` specifies the quadrant
that should be `m` x `n`.
"""
function part_2x2(A::Matrix, m::Int, n::Int, quad::String)
    if !(quad in ("TL", "TR", "BL", "BR"))
        throw(ArgumentError("""quad must be "TL", "TR", "BL", or "BR"."""))
    end

    hpart = (quad in ("TL", "TR")) ? m : size(A)[1] - m
    vpart = (quad in ("TL", "BL")) ? n : size(A)[2] - n

    TL, TR = A[1:hpart, 1:vpart], A[1:hpart, (vpart + 1):end]
    BL, BR = A[(hpart + 1):end, 1:vpart], A[(hpart + 1):end, (vpart + 1):end]

    return TL, TR, BL, BR
end

"""
    repart_1x2_to_1x3(AL::Matrix, AR::Matrix, n=1, side = "RIGHT")

Takes two submatrices and breaks the `side`hand portion into two
submatrices. `n` specifies the number of columns in what is ultimately
the middle submatrix.
"""
function repart_1x2_to_1x3(AL::Matrix, AR::Matrix, n=1, side = "RIGHT")
    if side == "RIGHT"
        vpart = n
        A1 = AL
        A2 = AR[:, 1:vpart]
        A3 = AR[:, (vpart + 1):end]
    else
        vpart = size(AL)[2] - n
        A1 = AL[:, 1:vpart]
        A2 = AL[:, (vpart + 1):end]
        A3 = AR
    end

    return A1, A2, A3
end

"""
    repart_2x1_to_3x1(xT::Vector, xB::Vector, m=1, side = "BOTTOM")

Repartition a vector already divided into top and bottom portions,
`xT` and `xB`, into three horizontal slabs, `x1`, `x2`, and `x3`.
The middle segment, `x2`, is created from the `m` elements in the `side`
input portion that are closest to the interface between `xT` and `xB`.

For example, if side == "BOTTOM", x1 == xT; x2 is created from the top m
elements from xB, and x2 is what remains of xB after its m elements are excluded.
"""
function repart_2x1_to_3x1(xT::Vector, xB::Vector, m=1, side = "BOTTOM")
    top_end = (side == "BOTTOM") ? size(xT)[1] : size(xT)[1] - m
    bottom_start = (side == "TOP") ? 1 : m + 1

    x1 = xT[1:top_end]
    x3 = xB[bottom_start: end]
    # What if you're taking from the bottom but the bottom is empty?
    # Or taking from the top, but the top is empty?
    # If "BOTTOM", x2 == xB[1:m] if length(xB) >= m,
    #              x2 == xB[1:0] if length(xB) == 0
    #              x2 == xB[1:length(xB)] if length(xB) <= m
    #              Generalization:
    #              x2 == xB[1:min(length(xB), m)]
    # What if m > length(xB) > 0?
    # If "TOP", x2 == xT[end - m + 1: end] if length(xT) >= m,
    #           x2 == xT[1:end] if length(xT) == 0
    #           x2 == xT[1:end] if length(xT) <= m
    #           Generalization:
    #           x2 == xT[max(1, end - m + 1):end]
    # What if m > length(xT) > 0?
    x2 = []
    if side == "BOTTOM"
        x2 = xB[1:min(bottom_start - 1, length(xB))]
    else
        x2 = xT[max(top_end + 1, 1): end]
    end
    # Expressed another way:
    # x2 = (side == "BOTTOM") ? xB[1:min(bottom_start - 1, length(xB))] : xT[max(top_end + 1, 1): end]
    return x1, x2, x3
end

"""
    repart_2x1_to_3x1(AT::Matrix, AB::Matrix, m=1, side = "BOTTOM")

Repartition a matrix already divided into top and bottom portions,
`AT` and `AB`, into three horizontal slabs, `A1`, `A2`, and `A3`.
The middle segment, `A2`, is created from the `m` rows in the `side`
input portion that are closest to the interface between `AT` and `AB`.

For example, if side == "BOTTOM", A1 == AT; A2 is created from the top m
rows from AB, and A3 is what remains of AB after its m rows are excluded.
"""
function repart_2x1_to_3x1(AT::Matrix, AB::Matrix, m=1, side = "BOTTOM")
    top_end = (side == "BOTTOM") ? size(AT)[1] : size(AT)[1] - m
    bottom_start = (side == "TOP") ? 1 : m + 1

    A1 = AT[1:top_end, :]
    A3 = AB[bottom_start: end, :]
    A2 = []
    if side == "BOTTOM"
        A2 = AB[1:min(bottom_start - 1, length(AB)), :]
    else
        A2 = AT[max(top_end + 1, 1): end, :]
    end
    # Expressed another way:
    # A2 = (side == "BOTTOM") ? AB[1:min(bottom_start - 1, length(AB)), :] : AT[max(top_end + 1, 1): end, :]
    return A1, A2, A3
end

"""
    repart_2x2_to_3x3(ATL::Matrix, ATR::Matrix,
                      ABL::Matrix, ABR::Matrix, m = 1, n = 1, quad = "BR")

Repartition a 2x2 matrix into a 3x3 matrix. The `quad` quadrant is broken into four
submatrices and the interior matrix returned, `A22`, will be m x n.
"""
function repart_2x2_to_3x3(ATL::Matrix, ATR::Matrix,
                           ABL::Matrix, ABR::Matrix, m = 1, n = 1, quad = "BR")
    hpart = (quad in ("TL", "TR")) ? size(ATL, 1) - m : m
    vpart = (quad in ("TL", "BL")) ? size(ATL, 2) - n : n

    if quad == "TL"
        A11, A12, A13 = ATL[1:hpart, 1:vpart], ATL[1:hpart, (vpart + 1):end], ATR[1:hpart, :]
        A21, A22, A23 = ATL[(hpart + 1):end, 1:vpart], ATL[(hpart + 1):end, (vpart + 1):end], ATR[(hpart + 1):end, :]
        A31, A32, A33 = ABL[:, 1:vpart], ABL[:, (vpart + 1):end], ABR[:, :]
    elseif quad == "TR"
        A11, A12, A13 = ATL[1:hpart, :], ATR[1:hpart, 1:vpart], ATR[1:hpart, (vpart + 1):end]
        A21, A22, A23 = ATL[(hpart + 1):end, :], ATR[(hpart + 1):end, 1:vpart], ATR[(hpart + 1):end, (vpart + 1):end]
        A31, A32, A33 = ABL[:, :], ABR[:, 1:vpart], ABR[:, (vpart + 1):end]
    elseif quad == "BL"
        A11, A12, A13 = ATL[:, 1:vpart], ATL[:, (vpart + 1):end], ATR[:, :]
        A21, A22, A23 = ABL[1:hpart, 1:vpart], ABL[1:hpart, (vpart + 1):end], ABR[1:hpart, :]
        A31, A32, A33 = ABL[(hpart + 1):end, 1:vpart], ABL[(hpart + 1):end, (vpart + 1):end], ABR[(hpart + 1):end, :]
    elseif quad == "BR"
        A11, A12, A13 = ATL[:, :], ATR[:, 1:vpart], ATR[:, (vpart + 1):end]
        A21, A22, A23 = ABL[1:hpart, :], ABR[1:hpart, 1:vpart], ABR[1:hpart, (vpart + 1):end]
        A31, A32, A33 = ABL[(hpart + 1):end, :], ABR[(hpart + 1):end, 1:vpart], ABR[(hpart + 1):end, (vpart + 1):end]
    end

    return ( A11, A12, A13,
             A21, A22, A23,
             A31, A32, A33 )

end

end
