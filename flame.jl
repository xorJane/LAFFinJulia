module flame

function part_2x1(A::Vector, ATsize=0, side="TOP")
    if ATsize < 0
        throw(DimensionMismatch("size < 0"))
    elseif ATsize > size(A)[1]
        throw(DimensionMismatch("size > row dimension"))
    elseif !(side in ("TOP", "BOTTOM"))
        throw(ArgumentError("""side must be "TOP" or "BOTTOM" """))
    end
    hpart = (side == "TOP") ? ATsize : size(A)[1] - ATsize
    AT, AB = A[1:hpart], A[hpart + 1:end]
end

function repart_2x1_to_3x1(AT::Vector, AB::Vector, m=1, side = "BOTTOM")
    top_end = (side == "BOTTOM") ? size(AT)[1] : size(AT)[1] - m
    bottom_start = (side == "TOP") ? 1 : m + 1

    A0 = AT[1:top_end]
    A2 = AB[bottom_start: end]
    # What if you're taking from the bottom but the bottom is empty?
    # Or taking from the top, but the top is empty?
    # If "BOTTOM", A1 == AB[1:m] if length(AB) >= m,
    #              A1 == AB[1:0] if length(AB) == 0
    #              A1 == AB[1:length(AB)] if length(AB) <= m
    #              Generalization:
    #              A1 == AB[1:min(length(AB), m)]
    # What if m > length(AB) > 0?
    # If "TOP", A1 == AT[end - m + 1: end] if length(AT) >= m,
    #           A1 == AT[1:end] if length(AT) == 0
    #           A1 == AT[1:end] if length(AT) <= m
    #           Generalization:
    #           A1 == AT[max(1, end - m + 1):end]
    # What if m > length(AT) > 0?
    A1 = []
    if side == "BOTTOM"
        A1 = AB[1:min(bottom_start - 1, length(AB))]
    else
        A1 = AT[max(top_end + 1, 1): end]
    end
    # Expressed another way:
    # A1 = (side == "BOTTOM") ? AB[1:min(bottom_start - 1, length(AB))] : AT[max(top_end + 1, 1): end]
    return A0, A1, A2
end

function cont_with_3x1_to_2x1(A0::Vector, A1::Vector, A2::Vector, side = "TOP")
    (side == "TOP") ? (vcat(A0, A1), A2) : (A0, vcat(A1, A2))
end

function merge_2x1!(T::Vector, B::Vector, A::Vector)
    Atemp = vcat(T, B)
    for i in 1:length(Atemp)
        A[i] = Atemp[i]
    end
end

end
