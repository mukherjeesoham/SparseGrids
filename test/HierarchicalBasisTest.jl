#-----------------------------------------------------
# SparseGrids | Hierarchical Basis Tests
# Soham 07/2020
#-----------------------------------------------------

lmax = 12
x = collect(range(0, stop=1, length=(2^lmax)+1))
y = NodalBasis(exp.(x))

@test HB_2_NB(NB_2_HB(y)).values ≈ y.values

ȳ = NB_2_HB(y)

using LinearAlgebra
function extract(HB::HierarchicalBasis{T})::Array{T,1} where {T}
    L2 = zeros(T, levels(HB))
    for level in 1:levels(HB)
        L2[level] = norm(HB.levels[level].values)
    end
    return L2
end

using PyPlot
plot(log10.(extract(ȳ)))
savefig("falloff.pdf")
close()


# 9 points --> exp(x) [Nodal Basis] --> Hierarchical basis (level 3)
# 100 points between [0, 1]
# y_analytic = exp(x) where x has 100 points
# y_numeric_upto_level_l = evaluate(HierarchichalBasis(up to level l)) on these 100 points
# error_upto_level_l = y_numeric_upto_level_l - y_analytic
# plot(upto_level_l, error_upto_level_l)

