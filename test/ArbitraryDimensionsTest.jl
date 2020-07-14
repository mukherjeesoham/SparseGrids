#-----------------------------------------------------
# SparseGrids | Arbitrary dimenions tests
# Soham 07/2020
#-----------------------------------------------------
# TODO: Test basis and x

D = 2
level = 2

U = map(NodalBasis{Float64,2}, 9, (x,y)->basis(2, CartesianIndex(4,4), (x,y)))
display(U.values)
