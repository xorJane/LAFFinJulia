include("../../flame.jl")
using .laff: trsv!

"""
    trsm_lnu!(L, B)

Solve LX = B, overwriting B with X, where L is a lower triangular matrix with a unit diagonal.
"""
function trsm_lnu!(L, B)

    BL, BR = flame.part_1x2(B, 
                            0, "LEFT")

    while size(BL, 2) < size(B, 2)

        B0, b1, B2 = flame.repart_1x2_to_1x3(BL, BR, 
                                             1, "RIGHT")

        #------------------------------------------------------------#

        laff.trsv!( "Lower triangular", "No transpose", "Unit diagonal", L, b1 )

        #------------------------------------------------------------#

        BL, BR = flame.cont_with_1x3_to_1x2(B0, b1, B2, 
                                            "LEFT")

    end

    flame.merge_1x2!(BL, BR, B)

end