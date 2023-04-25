abstract type Bravais_lattice end
struct Tetragonal <: Bravais_lattice end
struct Triangular <: Bravais_lattice end

include("triangular.jl")
