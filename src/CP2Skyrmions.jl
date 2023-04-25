module CP2Skyrmions
using LinearAlgebra, DelimitedFiles
using Random, NLopt, MPI, Printf

include("lattice_init.jl")
export Triangular, Jmat_init, neighbor_idx_table
include("cofig.jl")
export Paras
include("energy_fun.jl")
include("measurement.jl")
export getMagnetization, getEnergy
include("GradientDescent.jl")
export GDStatistics, GradientDescent

end
