module PracticeGemm

    include("PrintMMProblem.jl")

    function choose_problem_type()
        m = rand(2:4)
        k = rand(2:4)
        n = rand(2:4)
        tag = "Matrix Matrix Multiplication"
    
        case = rand(1:8)
        if case == 1
            m, k, n = 1, 1, 1
            tag = "Scalar Multiplication"
        elseif case == 2
            k, n = 1, 1
            tag = "SCAL"
        elseif case == 3
            m, k = 1, 1
            tag = "SCAL"
        elseif case == 4
            m, n = 1, 1
            tag = "DOT"
        elseif case == 5
            k = 1
            tag = "Outer Product"
        elseif case == 6
            n = 1
            tag = "Matrix-Vector Product"
        elseif case == 7
            m = 1
            tag = "Row Vector-Matrix Product"
        end
    
        return m, k, n, tag
    end

    function new_problem()
        m, k, n, tag = choose_problem_type()
        # Generate linear system
        A = rand(-2:2, m, k)
        B = rand(-2:2, k, n)
        # print linear system nicely
        println("$tag :")
        PrintMMProblem(A, B)
        return A, B
    end

    function show_answer(problem)
        A, B = problem
        C = A * B
        # print answer nicely
        PrintMMProblem(A, B, C)
    end

end
