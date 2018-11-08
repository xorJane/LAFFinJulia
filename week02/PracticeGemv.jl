module PracticeGemv

    include("./PrintMVProblem.jl")

    function new_problem()
        # Generate linear system
        m, n = rand(1:3, 2)
        A = rand(-2:2, m, n)
        x = rand(-2:2, n)
        # print linear system nicely
        PrintMVProblem(A, x)
        return A, x
    end

    function show_answer(problem)
        A, x = problem
        b = A * x
        # print answer nicely
        PrintMVProblem(A, x, b)
    end

end
