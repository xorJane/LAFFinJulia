include("./flame.jl")
using .flame
using Test

@testset "part_2x1 for Vector" begin
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

@testset "repart_2x1repart_2x1_to_3x1 for Vector" begin
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

@testset "cont_with_3x1_to_2x1 for Vectors" begin
    A0 = [1, 2, 3]
    A1 = [4]
    A2 = [5]
    @test flame.cont_with_3x1_to_2x1(A0, A1, A2) == ([1, 2, 3, 4], [5])
    @test flame.cont_with_3x1_to_2x1(A0, A1, A2, "BOTTOM") == ([1, 2, 3], [4, 5])
end

@testset "merge_2x1 for Vectors" begin
    yT = [1, 2, 3]
    yB = [4, 5]
    y = fill(0, 5)
    flame.merge_2x1!(yT, yB, y)
    @test y == [1, 2, 3, 4, 5]


end
