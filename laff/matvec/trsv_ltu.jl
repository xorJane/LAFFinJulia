using .laff: scal!, dots!

"""
    trsv_ltu!(L, b)

Solve L'x = b, overwriting b with x, where L is a lower triangular matrix with a unit diagonal.
"""
function trsv_ltu!(L, b)

    LTL, LTR, 
    LBL, LBR  = flame.part_2x2(L, 
                               0, 0, "TL")

    bT, 
    bB  = flame.part_2x1(b, 
                         0, "TOP")

    while size(LTL, 1) < size(L, 1)

        L00,  l01,      L02,  
        l10t, lambda11, l12t, 
        L20,  l21,      L22   = flame.repart_2x2_to_3x3(LTL, LTR, 
                                                        LBL, LBR, 
                                                        1, 1, "BR")

        b0,    
        beta1, 
        b2     = flame.repart_2x1_to_3x1(bT, 
                                         bB, 
                                         1, "BOTTOM")

        #------------------------------------------------------------#

        dots!( -l21, b2, beta1 )
#         scal!( 1.0 ./ lambda11, beta1 )

        #------------------------------------------------------------#

        LTL, LTR, 
        LBL, LBR  = flame.cont_with_3x3_to_2x2(L00,  l01,      L02,  
                                               l10t, lambda11, l12t, 
                                               L20,  l21,      L22,  
                                               "TL")

        bT, 
        bB  = flame.cont_with_3x1_to_2x1(b0,    
                                         beta1, 
                                         b2,    
                                         "TOP")

    end

    flame.merge_2x1!(bT, 
                     bB, b)

end