include("./laff.jl")
using Test
using LinearAlgebra

######### Test functions in util subdirectory ##########
@testset "util functions" begin
    # test onev with different input types!
    m = 3
    x1 = rand(m)
    laff.onev!(x1)
    @test x1 == ones(m)
    x2 = rand(m, 1); laff.onev!(x2)
    @test x2 == ones(m, 1)
    x3 = rand(1, m); laff.onev!(x3)
    @test x3 == ones(1, m)
    x1t = transpose(rand(m)); laff.onev!(x1t)
    @test x1t == transpose(ones(m))
    
    # test zerov with different input types!
    x1 = rand(m)
    laff.zerov!(x1)
    @test x1 == zeros(m)
    x2 = rand(m, 1); laff.zerov!(x2)
    @test x2 == zeros(m, 1)
    x3 = rand(1, m); laff.zerov!(x3)
    @test x3 == zeros(1, m)
    x1t = transpose(rand(m)); laff.zerov!(x1t)
    @test x1t == transpose(zeros(m))
end


######### Test functions in vecvec subdirectory ##########
@testset "vecvec functions" begin
    @testset "axpy!" begin
        ### Test axpy
        alpha1 = rand()
        alpha2 = [alpha1]
        m = 5
        x = rand(m)
        y = rand(m)
        # test axpy! method with numeric alpha
        updated_y1 = alpha1 * x + y
        laff.axpy!(alpha1, x, y)
        @test isapprox(y, updated_y1)

        # test axpy! method with array alpha
        updated_y2 = alpha2[1] * x + y
        laff.axpy!(alpha2, x, y)
        @test isapprox(y, updated_y2)
    end
    @testset "copy!" begin
        m = 5
        x, y = rand(m), rand(m)
        laff.copy!(x, y)
        @test isapprox(x, y)
    end
    @testset "dot and dots" begin
        # dot
        m = 5
        x, y = rand(m), rand(m)
        alpha = 5.0
        expected_alpha = alpha + x'y
        @test isapprox(laff.dot(x, y), x'y)
        @test isapprox(laff.dot(x, y, alpha), expected_alpha)
        
        alpha2 = [alpha]
        laff.dots!(x, y, alpha2)
        @test isapprox(alpha2[1], expected_alpha)
    end
    @testset "invscal and scal" begin
        # invscal
        m = 5
        alpha1 = rand(1, 1)
        alpha2 = alpha1[1]
        x1 = rand(m)
        x2 = copy(x1)
        y1 = copy(x1)
        y2 = copy(x1)
        expected_output = x1 ./ alpha1
        laff.invscal!(alpha1, x1)
        laff.invscal!(alpha2, x2)
        laff.scal!(1 ./alpha1, y1)
        laff.scal!(1 ./alpha2, y2)
        @test isapprox(expected_output, x1)
        @test isapprox(expected_output, x2)
        @test isapprox(expected_output, y1)
        @test isapprox(expected_output, y2)
    end
    @testset "norm2" begin
        x = rand(5)
        @test isapprox(laff.norm2(x), sqrt(sum(abs2.(x))))
    end
end

######### Test functions in matvec subdirectory ##########
@testset "matvec functions" begin
    @testset "gemv!" begin
        alpha = rand()
        m, n = 3, 4
        A = rand(m, n)
        x = rand(n)
        beta = rand(1, 1)
        y = rand(m)
        
        # test gemv! with no transpose
        trans = "No transpose"
        expected_output = alpha * A * x + beta[1, 1] * y
        laff.gemv!(trans, alpha, A, x, beta, y) 
        @test expected_output == y
        
        # test gemv! with no transpose
        trans = "Transpose"
        expected_output = beta[1,1] * A' * y + alpha * x
        laff.gemv!(trans, beta, A, y, alpha, x)
        @test expected_output == x
        # Note: above, we also made sure that alpha and beta can both be either
        # scalars or Matrices containing a single scalar
    end
    
    @testset "ger!" begin
        m, n = 3, 4
        x = rand(n)
        y = rand(m)
        
        # test ger! where alpha is a scalar
        alpha = rand()
        A = rand(m, n)
        expected_output = alpha * y * x' + A
        laff.ger!(alpha, y, x, A)
        @test A == expected_output
        
        # test ger! where alpha is a Matrix
        alpha = rand(1, 1)
        A = rand(m, n)
        expected_output = alpha[1,1] * y * x' + A
        laff.ger!(alpha, y, x, A)
        @test A == expected_output
    end
    
    @testset "trsv*! functions" begin
        # test trsv_lnn! alone
        m = 4
        L = [i >= j ? rand() : 0.0 for i in 1:m, j in 1:m]
        x = rand(m)
        b = L * x
        laff.trsv_lnn!(L, b)
        @test isapprox(x, b)
        
        # test trsv_lnn! through trsv!
        uplo = "Lower triangular"
        trans = "No transpose"
        diag = "Nonunit diagonal"
        x = rand(m)
        b = L * x
        laff.trsv!(uplo, trans, diag, L, b)
        @test isapprox(x, b)
        ######################
        # test trsv_lnu! alone
        L = [i > j ? rand() : ((i == j) ? 1.0 : 0.0) for i in 1:m, j in 1:m]
        x = rand(m)
        b = L * x
        laff.trsv_lnu!(L, b)
        @test isapprox(x, b)
        
        # test trsv_lnu! through trsv!
        uplo = "Lower triangular"
        trans = "No transpose"
        diag = "Unit diagonal"
        x = rand(m)
        b = L * x
        laff.trsv!(uplo, trans, diag, L, b)
        @test isapprox(x, b)
        
        ######################
        # test trsv_unn! alone
        U = [i <= j ? rand() : 0.0 for i in 1:m, j in 1:m]
        x = rand(m)
        b = U * x
        laff.trsv_unn!(U, b)
        @test isapprox(x, b)
        
        # test trsv_unn! through trsv!
        uplo = "Upper triangular"
        trans = "No transpose"
        diag = "Nonunit diagonal"
        A = rand(m, m)
        U = UpperTriangular(A)
        x = rand(m)
        b = U * x
        laff.trsv!(uplo, trans, diag, A, b)
        @test isapprox(x, b)
        
        # test trsv_unu! alone
        # Give `A` a unit diagonal
        A = [i==j ? 1.0 : A[i, j] for i in 1:m, j in 1:m]
        U = UpperTriangular(A)
        b = U * x
        laff.trsv_unu!(A, b)
        @test isapprox(x, b)
        
        # test trsv_unu! through trsv!
        uplo = "Upper triangular"
        trans = "No transpose"
        diag = "Unit diagonal"
        b = U * x
        laff.trsv!(uplo, trans, diag, A, b)
        @test isapprox(x, b)
        
        # test trsv_utn!
        m = 4
        A = rand(m, m)
        U = UpperTriangular(A)
        x = rand(m)
        
        # test trsv_utn! alone
        b = U' * x
        laff.trsv_utn!(A, b)
        @test isapprox(x, b)
        
        # test trsv_utn! through trsv!
        uplo = "Upper triangular"
        trans = "Transpose"
        diag = "Nonunit diagonal"
        b = U' * x
        laff.trsv!(uplo, trans, diag, A, b)
        @test isapprox(x, b)   
        
        ######################
        # test trsv_ltu! alone
        
        # Assume function works by FIRST extracting the
        # lower triangular part  and then transposing A
        for i in 1:m; A[i, i] = 1.0; end
        L = LowerTriangular(A)
        b = L' * x
        laff.trsv_ltu!(A, b)
        @test isapprox(x, b)
        
        # test trsv_ltu! through trsv!
        uplo = "Lower triangular"
        trans = "Transpose"
        diag = "Unit diagonal"
        b = L' * x
        laff.trsv!(uplo, trans, diag, A, b)
        @test isapprox(x, b) 
    end
end 

######### Test functions in matmat subdirectory ##########
@testset "matmat functions" begin
    # Test gemm
    m, n, k = 4, 5, 6
    C = rand(m, n)
    A = rand(m, k)
    B = rand(k, n)
    alpha = rand()
    beta = rand(1)
    expected_C = beta[1] * C + alpha * A * B
    laff.gemm!(alpha, A, B, beta, C)
    @test isapprox(C, expected_C)
     
    @testset "trsm!" begin
        m = 4
        A = rand(m, m)
#         A = [i + m*(j - 1) for i in 1:m, j in 1:m]
        U = UpperTriangular(A)
#         X = A .+ ones(m, m)
        X = rand(m, m)
        #################
        # Test trsm_unn!
        # Test trsm_unn! alone
        B = U * X
        laff.trsm_unn!(A, B)
        @test isapprox(B, X)
        
        # Test trsm_unn! from trsm!
        uplo = "Upper triangular"
        trans = "No transpose"
        diag = "Nonunit diagonal"
        B = U * X
        laff.trsm!(uplo, trans, diag, A, B)
        @test isapprox(B, X)
        
        
        #################
        # Test trsm_lnu!
        # Test trsm_lnu! alone
        for i in 1:size(A, 1); A[i, i] = 1.0; end
        L = LowerTriangular(A)
        B = L * X
        laff.trsm_lnu!(A, B)
        @test isapprox(B, X)

        # Test trsm_lnu! from trsm!
        uplo = "Lower triangular"
        trans = "No transpose"
        diag = "Unit diagonal"
        B = L * X
        laff.trsm!(uplo, trans, diag, A, B)
        @test isapprox(B, X)
        
        m = 4
        A = rand(m, m)
        U = UpperTriangular(A)
        X = rand(m, m)
        #################
        # Test trsm_utn!
        # Test trsm_utn! alone
        B = U' * X
        laff.trsm_utn!(A, B)
        @test isapprox(B, X)
        
        # Test trsm_utn! from trsm!
        uplo = "Upper triangular"
        trans = "Transpose"
        diag = "Nonunit diagonal"
        B = U' * X
        laff.trsm!(uplo, trans, diag, A, B)
        @test isapprox(B, X)
        
        # Test trsm_ltu!
        # Test trsm_ltu! alone
        for i in 1:m; A[i, i] = 1.0; end
        L = LowerTriangular(A)
        B = L' * X
        laff.trsm_ltu!(A, B)
        @test isapprox(B, X)
        
        # Test trsm_ltu! from trsm!
        uplo = "Lower triangular"
        trans = "Transpose"
        diag = "Unit diagonal"
        B = L' * X
        laff.trsm!(uplo, trans, diag, A, B)
        @test isapprox(B, X) 
     
    end
end