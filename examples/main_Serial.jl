using CP2Skyrmions

D_ion = 15.0  # single-ion anisotropy
H_min = 8.0   # min value of h
H_max = 12.0  # max value of h
H_num = 20    # number of magnetic field

# construct the J1-J2 triangular lattice
J_intra = [-1.0, 1.0, 2.6, 0.0, 0.0, 0.0,
2/(1+√5), 1.0, 2.6, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# system size
N1, N2 = 10, 10
Ntot = N1 * N2

max_dist = 3
lat = Triangular()
J_mat = Jmat_init(J_intra, max_dist, lat)
basis = 1
Nei_table = neighbor_idx_table(N1, N2, max_dist, lat, basis)
N_bond = 18

# the vector of magnetic fields
Hz = range(H_min, stop=H_max, length=H_num)

for h_ext in Hz
	gll = Paras(J_mat, D_ion, h_ext, N1, N2, Ntot, N_bond, Nei_table)
	datadir = joinpath(@__DIR__, "data", string(N1)*"_"*string(N2), "D_"*string(D_ion), "H_"*string(round(h_ext, digits=3)))
	# number of random seeds per calculation
	numRand = 50
	# number of iterations
	numIter = 1
	gds = GDStatistics(numRand, 0, 1, 0, numIter, 0, 0)
	GradientDescent(gll, gds, datadir)
end