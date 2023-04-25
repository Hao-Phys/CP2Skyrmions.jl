using CP2Skyrmions, MPI, DelimitedFiles

@assert length(ARGS) >= 3
const D_ion = parse(Float64, ARGS[1])
const H_min = parse(Float64, ARGS[2])
const H_max = parse(Float64, ARGS[3])
length(ARGS) == 4 && const D_prev = ARGS[4]


function main_MPI()
	MPI.Initialized() || MPI.Init()
	commSize = MPI.Comm_size(MPI.COMM_WORLD)
	commRank = MPI.Comm_rank(MPI.COMM_WORLD)

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

	Hz = range(H_min, stop=H_max, length=commSize)[commRank+1]
	gll = Paras(J_mat, D_ion, Hz, N1, N2, Ntot, N_bond, Nei_table)
	numRand, numIter = 50, 2
	gds = GDStatistics(numRand, 0, commSize, 0, numIter, 0, 0)
	datadir = joinpath(@__DIR__, "data/J1-J2", string(N1)*"_"*string(N2), "D_"*string(D_ion), "H_"*string(round(Hz, digits=3)))

	if length(ARGS) == 4
		fname = joinpath(@__DIR__, "data/J1-J2", string(N1)*"_"*string(N2), "D_"*string(D_prev), "opt_angles.dat")
		if isfile(fname)
			println("use previous data of D= ", D_prev)
			prevData = readdlm(fname)
			gds.numPrev = length(prevData[1, :])
			GradientDescent(gll, gds, datadir, prevData)
		else
			println("previous data of D_", D_prev, "/opt_angles.dat does not exist! Run the simulation without previous data")
			GradientDescent(gll, gds, datadir)
		end
	else
		GradientDescent(gll, gds, datadir)
	end
end

@time main_MPI()
