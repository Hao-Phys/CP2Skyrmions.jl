mutable struct Paras
	J_mat :: Array{Float64, 3}
	D_ion :: Float64
	h_ext :: Float64
	# lattice size
	N1   :: Int
	N2   :: Int
	Ntot :: Int
	# number of exchange bonds
	N_bond :: Int
	# nearest-neighbors table
	Nei_table :: Array{Int, 2}
end

function coor2site(gll::Paras, a, b)
	idx = a + (b-1) * gll.N1
end

function site2coor(gll::Paras, idx)
	b = fld(idx-1, gll.N1) + 1
	a = mod(idx-1, gll.N1) + 1
	a, b
end
