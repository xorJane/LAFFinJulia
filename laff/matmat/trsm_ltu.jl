using .laff: trsv!

"""
    trsm_ltu!(L, B)

Solve L'X = B, overwriting B with X, where L is a lower triangular matrix with a unit diagonal.
"""
function trsm_ltu!(L, B)

    BT, 
    BB  = flame.part_2x1(B, 
                         0, "TOP")

    while size(BT, 1) < size(B, 1)

        B0,  
        b1t, 
        B2   = flame.repart_2x1_to_3x1(BT, 
                                       BB, 
                                       1, "BOTTOM")

        #------------------------------------------------------------#

        laff.trsv!( "Lower triangular", "Transpose", "Unit diagonal", L, b1t )

        #------------------------------------------------------------#

        BT, 
        BB  = flame.cont_with_3x1_to_2x1(B0,  
                                         b1t, 
                                         B2,  
                                         "TOP")

    end

    flame.merge_2x1!(BT, 
                     BB, B)

end