include("./flame.jl")
using .flame
using Test

@testset "merge_2x2" begin
    # Divide square matrix into four evenly
    # sized quadrants and reconstruct
    A = reshape([x for x in 1:16], 4, 4)
    Anew = fill(0, 4, 4)
    TL = A[1:2, 1:2]
    TR = A[1:2, 3:4]
    BL = A[3:4, 1:2]
    BR = A[3:4, 3:4]
    flame.merge_2x2!(TL, TR, BL, BR, Anew) 
    @test Anew == A
    ## Automate testing for matrices of various
    #  sizes
    for m in 2:5, n in 2:5
        A = reshape([x for x in 1:m*n], m, n)
        Anew = fill(0, m, n)
        # pick horizontal and vertical lines
        # along which to divide A into quadrants.
        for hsplit in 1:(m - 1), vsplit in 1:(n - 1)
            TL = A[1:hsplit, 1:vsplit]
            TR = A[1:hsplit, (vsplit + 1):end]
            BL = A[(hsplit + 1):end, 1:vsplit]
            BR = A[(hsplit + 1):end, (vsplit + 1):end]
            flame.merge_2x2!(TL, TR, BL, BR, Anew) 
            @test Anew == A
        end
    end
end

@testset "merge_2x1" begin
    @testset "for Vectors" begin
        yT = [1, 2, 3]
        yB = [4, 5]
        y = fill(0, 5)
        flame.merge_2x1!(yT, yB, y)
        @test y == [1, 2, 3, 4, 5]
    end

    @testset "for Matrices" begin
        AT = [1 2 3;]
        AB = [4 5 6; 7 8 9]
        A = fill(0, 3, 3)
        flame.merge_2x1!(AT, AB, A)
        @test A == [i + 3*(j - 1) for j in 1:3, i in 1:3]
    end
end

@testset "merge_1x2!" begin
    A = reshape([x for x in 1:25], 5, 5)
    Anew = fill(0, 5, 5)
    L = A[:, 1:2]
    R = A[:, 3:5]
    flame.merge_1x2!(L, R, Anew)
    @test A == Anew
    Anew = fill(0, 1, 5)
    L = A[1:1, 1:2]
    R = A[1:1, 3:5]
    flame.merge_1x2!(L, R, Anew)
    @test A[1:1, :] == Anew
end

@testset "cont_with_1x3_to_1x2" begin
    @testset "for Matrices" begin
        m, n = 3, 4
        A = [(j + n*(i - 1)) for i in 1:m, j in 1:n]
        A0 = A[:, 1:2]
        A1 = A[:, 3:3]
        A2 = A[:, 4:end]
        @test flame.cont_with_1x3_to_1x2(A0, A1, A2) == ([1 2 3; 5 6 7; 9 10 11], A2)
        @test flame.cont_with_1x3_to_1x2(A0, A1, A2, "RIGHT") == (A0, [3 4; 7 8; 11 12])
    end
end

@testset "cont_with_3x1_to_2x1" begin
    @testset "for Vectors" begin
        A0 = [1, 2, 3]
        A1 = [4]
        A2 = [5]
        @test flame.cont_with_3x1_to_2x1(A0, A1, A2) == ([1, 2, 3, 4], [5])
        @test flame.cont_with_3x1_to_2x1(A0, A1, A2, "BOTTOM") == ([1, 2, 3], [4, 5])
    end
    @testset "for Matrices" begin
        A0 = [1 2 3;]
        A1 = [4 5 6;]
        A2 = [7 8 9;]
        @test flame.cont_with_3x1_to_2x1(A0, A1, A2) == ([1 2 3; 4 5 6], A2)
        @test flame.cont_with_3x1_to_2x1(A0, A1, A2, "BOTTOM") == (A0, [4 5 6; 7 8 9])
    end
end

@testset "cont_with_3x3_to_2x2" begin
    A = [1 2 3; 4 5 6; 7 8 9]
    A00, A01, A02 = A[1:1, 1:1], A[1:1, 2:2], A[1:1, 3:3]
    A10, A11, A12 = A[2:2, 1:1], A[2:2, 2:2], A[2:2, 3:3]
    A20, A21, A22 = A[3:3, 1:1], A[3:3, 2:2], A[3:3, 3:3]
    # Test each configuration of recombining
    TL, TR, BL, BR = flame.cont_with_3x3_to_2x2(A00, A01, A02, A10, A11, A12, A20, A21, A22,"TL")
    @test (TL, TR, BL, BR) == (A[1:2, 1:2], A[1:2, 3:3], A[3:3, 1:2], A[3:3, 3:3])
    TL, TR, BL, BR = flame.cont_with_3x3_to_2x2(A00, A01, A02, A10, A11, A12, A20, A21, A22,"TR")
    @test (TL, TR, BL, BR) == (A[1:2, 1:1], A[1:2, 2:3], A[3:3, 1:1], A[3:3, 2:3]) 
    TL, TR, BL, BR = flame.cont_with_3x3_to_2x2(A00, A01, A02, A10, A11, A12, A20, A21, A22,"BL")
    @test (TL, TR, BL, BR) == (A[1:1, 1:2], A[1:1, 3:3], A[2:3, 1:2], A[2:3, 3:3])
    TL, TR, BL, BR = flame.cont_with_3x3_to_2x2(A00, A01, A02, A10, A11, A12, A20, A21, A22,"BR")
    @test (TL, TR, BL, BR) == (A[1:1, 1:1], A[1:1, 2:3], A[2:3, 1:1], A[2:3, 2:3])
    
end

@testset "part 1x2 for Matrix" begin
    A = [1 2 3; 4 5 6; 7 8 9]
    @test_throws DimensionMismatch flame.part_1x2(A, -1)
    @test_throws DimensionMismatch flame.part_1x2(A, size(A)[2] + 1)
    @test_throws ArgumentError flame.part_1x2(A, 0, "HELLO")
    @test flame.part_1x2(A, 0, "LEFT") == (Array{Int64}(undef, 3, 0), A)
    @test flame.part_1x2(A, 0, "RIGHT") == (A, Array{Int64}(undef, 3, 0))
    @test flame.part_1x2(A, 1, "LEFT") == (reshape(A[:, 1], 3, 1), A[:, 2:3])
    @test flame.part_1x2(A, 1, "RIGHT") == (A[:, 1:2], reshape(A[:, 3], 3, 1))
    @test flame.part_1x2(A, 3, "LEFT") == (A, Array{Int64}(undef, 3, 0))
    @test flame.part_1x2(A, 3, "RIGHT") == (Array{Int64}(undef, 3, 0), A )
end

@testset "part_2x1" begin
    @testset "for Vectors" begin
        x = [1, 2, 3, 4]
        @test_throws DimensionMismatch flame.part_2x1(x, -1)
        @test_throws DimensionMismatch flame.part_2x1(x, 5)
        @test_throws ArgumentError flame.part_2x1(x, 0, "HELLO")
        @test flame.part_2x1(x, 0, "TOP") == ([], x)
        @test flame.part_2x1(x, 0, "BOTTOM") == (x, [])
        @test flame.part_2x1(x, 1, "TOP") == (x[1:1], x[2:end])
        @test flame.part_2x1(x, 1, "BOTTOM") == (x[1:end - 1], x[end:end])
        @test flame.part_2x1(x, 2, "TOP") == (x[1:2], x[3:end])
        @test flame.part_2x1(x, 2, "BOTTOM") == (x[1:end - 2], x[end - 1:end])
        @test flame.part_2x1(x, 4, "TOP") == (x[1:end], [])
        @test flame.part_2x1(x, 4, "BOTTOM") == ([], x[1:end])
    end
    @testset "for Matrices" begin
        A = reshape([x for x in 1:16], 4, 4)
        @test_throws DimensionMismatch flame.part_2x1(A, -1)
        @test_throws DimensionMismatch flame.part_2x1(A, 5)
        @test_throws ArgumentError flame.part_2x1(A, 0, "HELLO")
        @test flame.part_2x1(A, 0, "TOP") == (Array{Int64}(undef, 0, 4), A)
        @test flame.part_2x1(A, 0, "BOTTOM") == (A, Array{Int64}(undef, 0, 4),)
        @test flame.part_2x1(A, 1, "TOP") == (A[1:1, :], A[2:end, :])
        @test flame.part_2x1(A, 1, "BOTTOM") == (A[1:end - 1, :], A[end:end, :])
        @test flame.part_2x1(A, 2, "TOP") == (A[1:2, :], A[3:end, :])
        @test flame.part_2x1(A, 2, "BOTTOM") == (A[1:end - 2, :], A[end - 1:end, :])
        @test flame.part_2x1(A, 4, "TOP") == (A[1:end, :], Array{Int64}(undef, 0, 4))
        @test flame.part_2x1(A, 4, "BOTTOM") == (Array{Int64}(undef, 0, 4), A[1:end, :])
    end
end

@testset "part_2x2" begin
    m, n = 5,5
    A = reshape([x for x in 1:m*n], m, n)
    for m in 0:(m - 1), n in 0:(n - 1)
        TL, TR, BL, BR = flame.part_2x2(A, m, n, "TL")
        @test size(TL) == (m, n)
        @test size(BR) == size(A) .- (m, n)
        @test size(TR) == (m, size(A)[2] - n)
        @test size(BL) == (size(A)[1] - m, n)
    end
    
    for m in 0:(m - 1), n in 0:(n - 1)
        TL, TR, BL, BR = flame.part_2x2(A, m, n, "TR")
        @test size(TL) == (m, size(A)[2] - n)
        @test size(BR) == (size(A)[1] - m, n)
        @test size(TR) == (m, n)
        @test size(BL) == size(A) .- (m, n)
    end
    
    for m in 0:(m - 1), n in 0:(n - 1)
        TL, TR, BL, BR = flame.part_2x2(A, m, n, "BL")
        @test size(TL) == (size(A)[1] - m, n)
        @test size(BR) == (m, size(A)[2] - n)
        @test size(TR) == size(A) .- (m, n)
        @test size(BL) == (m, n)
    end
    
    for m in 0:(m - 1), n in 0:(n - 1)
        TL, TR, BL, BR = flame.part_2x2(A, m, n, "BR")
        @test size(TL) == size(A) .- (m, n)
        @test size(BR) == (m, n)
        @test size(TR) == (size(A)[1] - m, n)
        @test size(BL) == (m, size(A)[2] - n)
    end
end

@testset "repart_1x2_to_1x3" begin
    A = [1 2 3; 4 5 6; 7 8 9]
    n = 1
    @test flame.repart_1x2_to_1x3(A[:, 1:1], A[:, 2:3], n, "RIGHT") == (A[:, 1:1], A[:, 2:2], A[:, 3:3])
    @test flame.repart_1x2_to_1x3(A[:, 1:2], A[:, 3:3], n, "LEFT") == (A[:, 1:1], A[:, 2:2], A[:, 3:3])
    # What if one of the portions of A is too small to further subdivide?
    @test flame.repart_1x2_to_1x3(A[:, 1:1], A[:, 2:3], n, "LEFT") == (Array{Int64}(undef, 3, 0), A[:, 1:1], A[:, 2:3])
    @test flame.repart_1x2_to_1x3(A[:, 1:2], A[:, 3:3], n, "RIGHT") == (A[:, 1:2], A[:, 3:3], Array{Int64}(undef, 3, 0))
end

@testset "repart_2x1_to_3x1" begin
    @testset "for Vectors" begin
        x = [1, 2, 3, 4]
        for ATsize in 0:length(x)
            AT = x[1:ATsize]
            AB = x[ATsize + 1:end]

            m = 1
            # For "bottom"
            A0 = AT[1:ATsize]
            A1 = AB[1:min(m, length(AB))]
            A2 = AB[m + 1:end]
            @test flame.repart_2x1_to_3x1(AT, AB, m, "BOTTOM") == (A0, A1, A2)
            # For "top"
            A0 = AT[1:ATsize - m]
            A1 = AT[max(ATsize - m + 1, 1):end]
            A2 = AB[1:end]
            @test flame.repart_2x1_to_3x1(AT, AB, m, "TOP") == (A0, A1, A2)
        end
    end
    @testset "for Matrices" begin
        A = reshape([i for i in 1:16], 4, 4)
        for ATsize in 0:size(A)[1]
            AT = A[1:ATsize, :]
            AB = A[ATsize + 1:end, :]

            m = 1
            # For "bottom"
            A0 = AT[1:ATsize, :]
            A1 = AB[1:min(m, length(AB)), :]
            A2 = AB[m + 1:end, :]
            @test flame.repart_2x1_to_3x1(AT, AB, m, "BOTTOM") == (A0, A1, A2)
            # For "top"
            A0 = AT[1:ATsize - m, :]
            A1 = AT[max(ATsize - m + 1, 1):end, :]
            A2 = AB[1:end, :]
            @test flame.repart_2x1_to_3x1(AT, AB, m, "TOP") == (A0, A1, A2)
        end
    end
end

@testset "repart_2x2_to_3x3" begin
    A = reshape([i for i in 1:16], 4, 4)
    ATL, ATR = A[1:2, 1:2], A[1:2, 3:4]
    ABL, ABR = A[3:4, 1:2], A[3:4, 3:4]
    output_matrices = flame.repart_2x2_to_3x3(ATL, ATR, ABL, ABR, 1, 1, "TL")
    @test output_matrices[1:3] == (ATL[1:1, 1:1], ATL[1:1, 2:2], ATR[1:1, :])
    @test output_matrices[4:6] == (ATL[2:2, 1:1], ATL[2:2, 2:2], ATR[2:2, :])
    @test output_matrices[7:9] == (ABL[:, 1:1], ABL[:, 2:2], ABR[:, :])
    output_matrices = flame.repart_2x2_to_3x3(ATL, ATR, ABL, ABR, 1, 1, "TR")
    @test output_matrices[1:3] == (ATL[1:1, :], ATR[1:1, 1:1], ATR[1:1, 2:2])
    @test output_matrices[4:6] == (ATL[2:2, :], ATR[2:2, 1:1], ATR[2:2, 2:2])
    @test output_matrices[7:9] == (ABL[:, :], ABR[:, 1:1], ABR[:, 2:2])
    output_matrices = flame.repart_2x2_to_3x3(ATL, ATR, ABL, ABR, 1, 1, "BL")
    @test output_matrices[1:3] == (ATL[:, 1:1], ATL[:, 2:2], ATR[:, :])
    @test output_matrices[4:6] == (ABL[1:1, 1:1], ABL[1:1, 2:2], ABR[1:1, :])
    @test output_matrices[7:9] == (ABL[2:2, 1:1], ABL[2:2, 2:2], ABR[2:2, :])
    output_matrices = flame.repart_2x2_to_3x3(ATL, ATR, ABL, ABR, 1, 1, "BR")
    @test output_matrices[1:3] == (ATL[:, :], ATR[:, 1:1], ATR[:, 2:2])
    @test output_matrices[4:6] == (ABL[1:1, :], ABR[1:1, 1:1], ABR[1:1, 2:2])
    @test output_matrices[7:9] == (ABL[2:2, :], ABR[2:2, 1:1], ABR[2:2, 2:2])
end

