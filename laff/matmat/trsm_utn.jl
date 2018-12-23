using .laff: trsv!

"""
    trsm_utn!(U, B)

Solve U'X = B, overwriting B with X, where U is an upper triangular matrix with a nonunit diagonal.
"""
function trsm_utn!(U, B)

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

        laff.trsv!( "Upper triangular", "Transpose", "Nonunit diagonal", U, b1t )

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