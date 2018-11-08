using Printf

function PrintMVProblem( A, x, b = false )
    let
        m_A, n_A = size( A )
        m_x = length(x)
        middlerow = max(1, max(cld(m_A, 2), cld(m_x, 2)))

        for i = 1:max( m_A, m_x )
            outstring = ""
            # Leftmost character -- create matrix bracket
            if m_A == 1 && i == 1
                outstring *= "< "
            elseif i == 1
                outstring *= "/ "
            elseif i == m_A
                outstring *= "\\ "
            elseif i < m_A
                outstring *= "| "
            else
                outstring *= "  "
            end

            # Add matrix elements for ith row to output string
            for j = 1:n_A
                if i <= m_A
                    outstring *= @sprintf( "%2d ", A[ i,j ] )
                else
                    outstring *= "   "
                end
            end

            # Add right matrix bracket
            if m_A == 1 && i == 1
                outstring *= " >  "
            elseif i == 1
                outstring *= " \\  "
            elseif i == m_A
                outstring *= " /  "
            elseif i < m_A
                outstring *= " |  "
            else
                outstring *= "    "
            end

            # Add left bracket for Vector
            if m_x == 1 && i== 1
                outstring *= "< "
            elseif i == 1
                outstring *= "/ "
            elseif i == m_x
                outstring *= "\\ "
            elseif i < m_x
                outstring *= "| "
            else
                outstring *= "  "
            end

            # print element in the current row from vector x
            if i <= m_x
                outstring *= @sprintf( "%2d ", x[i] )
            else
                outstring *= "   "
            end

            # Add right bracket for vector x
            if m_x == 1 && i == 1
                outstring *= ">  "
            elseif i == 1
                outstring *= "\\  "
            elseif i == m_x
                outstring *= "/  "
            elseif i < m_x
                outstring *= "|  "
            else
                outstring *= "   "
            end

            # Add equal sign to equation
            if i == middlerow
                outstring *= " =  "
            else
                outstring *= "    "
            end

            # if b is passed (optional input),
            # add it to the outstring as well
            if b != false
                m_b = length(b)

                # Add left bracket for Vector, b
                if m_b == 1 && i== 1
                    outstring *= "< "
                elseif i == 1
                    outstring *= "/ "
                elseif i == m_b
                    outstring *= "\\ "
                elseif i < m_b
                    outstring *= "| "
                else
                    outstring *= "  "
                end

                # print element in the current row from vector b
                if i <= m_b
                    outstring *= @sprintf( "%2d ", b[i] )
                end

                # Add right bracket for vector b
                if m_b == 1 && i == 1
                    outstring *= ">  "
                elseif i == 1
                    outstring *= "\\  "
                elseif i == m_b
                    outstring *= "/  "
                elseif i < m_b
                    outstring *= "|  "
                else
                    outstring *= "   "
                end
            end

            println(outstring)
        end
    end
end
