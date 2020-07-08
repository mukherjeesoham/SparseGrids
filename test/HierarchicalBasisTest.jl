#-----------------------------------------------------
# SparseGrids | Hierarchical Basis Tests
# Soham 07/2020
#-----------------------------------------------------

lmax = 2
x = collect(range(0, stop=1, length=(2^lmax)+1))
y = NodalBasis(exp.(x))

@test HB_2_NB(NB_2_HB(y)).values ≈ y.values
