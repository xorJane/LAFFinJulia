using .laff: trsv!

"""
    trsm_utn!(U, B)

Solve U'X = B, overwriting B with X, where U is an upper triangular matrix with a nonunit diagonal.
"""
function trsm_utn!(U, B)

    BL, BR = flame.part_1x2(B, 
                            0, "LEFT")

    while size(BL, 2) < size(B, 2)

        B0, b1, B2 = flame.repart_1x2_to_1x3(BL, BR, 
                                             1, "RIGHT")

        #------------------------------------------------------------#

        laff.trsv!( "Upper triangular", "Transpose", "Nonunit diagonal", U, b1 )

        #------------------------------------------------------------#

        BL, BR = flame.cont_with_1x3_to_1x2(B0, b1, B2, 
                                            "LEFT")

    end

    flame.merge_1x2!(BL, BR, B)

end