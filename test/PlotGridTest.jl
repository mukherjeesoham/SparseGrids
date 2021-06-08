#-----------------------------------------------------
# SparseGrids | Plot grid points 
# Soham 07/2020
#-----------------------------------------------------

using PyPlot

l = 2 
N = 2^l + 1 

for i in 1:N, j in 1:N
    x = i
    y = j
    
    # Nodal
    plot(x, y, "k o", fillstyle="none", markersize=10)

    # level 0, 0
    if (i == 1 || i == N) && (j == 1 || j == N)
        plot(x, y, "b o")
    end

        # level 1, 0
        if (i == 3)  && (j == 1 || j == N)
            plot(x, y, "y o")
        end

        # level 0, 1
        if (j == 3)  && (i == 1 || i == N)
            plot(x, y, "y o")
        end
        
        # level 1, 1
        if (j == 3)  && (i == 3)
            plot(x, y, "r o")
        end

        # level 2, 0
        if (i == 2 || i == 4)  && (j == 1 || j == N)
            plot(x, y, "r o")
        end

    if false
        # level 0, 2
        if (j == 2 || j == 4)  && (i == 1 || i == N)
            plot(x, y, "r o")
        end

        # level 2, 1
        if (i == 2 || i == 4)  && (j == 3)
            plot(x, y, "c o")
        end

        # level 1, 2
        if (j == 2 || j == 4)  && (i == 3)
            plot(x, y, "c o")
        end

        # level 2, 2
        if (j == 2 || j == 4)  && (i == 2 || i == 4)
            plot(x, y, "m o")
        end
    end
end

xlim(0.7, 5.3)
ylim(0.7, 5.3)
savefig("sequence/img_13.pdf")
close()
    

