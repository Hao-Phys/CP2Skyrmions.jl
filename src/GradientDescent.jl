mutable struct GDStatistics
	numRand :: Int
	numRacc :: Int
	numGlob :: Int
	numGacc :: Int
	numIter :: Int
	numPrev :: Int
	numPacc :: Int
	# GDStatistics() = new(0, 0, 0, 0, 0)
end

"""
write to file from time to time, in case the program crashes, one can still get
some intermediate results.
"""
function WriteResults(
	gll :: Paras, gds :: GDStatistics,
	opt_angles :: Vector{Float64},
	datadir :: String, Iter :: Int
	)

	O_vals = zeros(Float64, 4*gll.Ntot)
	O_site = zeros(Float64, 4)

	for ii = 1:gll.Ntot
		θ  = opt_angles[(ii-1)*4+1]
		φ  = opt_angles[(ii-1)*4+2]
		α₁ = opt_angles[(ii-1)*4+3]
		α₂ = opt_angles[(ii-1)*4+4]
		O_ave4!(θ, φ, α₁, α₂, O_site)
		O_vals[(ii-1)*4+1:(ii-1)*4+4] = O_site
	end

	fname1 = joinpath(datadir, "angles_run"*string(Iter)*".dat")
	writedlm(fname1, opt_angles)
	fname2 = joinpath(datadir, "O_vals_run"*string(Iter)*".dat")
	writedlm(fname2, O_vals)

	opt_energy = getEnergy(O_vals, gll)
	opt_magnet = getMagnetization(O_vals, gll)
	summary = [gll.h_ext, opt_energy, opt_magnet]
	fname3 = joinpath(datadir, "summary_run"*string(Iter)*".dat")
	writedlm(fname3, summary)
	stats = [gds.numRand, gds.numRacc, gds.numGlob, gds.numGacc]
	fname4 = joinpath(datadir, "stat_run"*string(Iter)*".dat")
	writedlm(fname4, stats)
end

"""
The gradient descent function.
`gll`: object (struct) defined in cofig.jl
`gds`: object (struct) for recording the statistics of the calculation
`datadir`: the path to save data
`prevData` (optional argument): previous data for improving the quality of minimization, default to be empty
`saveRate` (optional argument): the rate to write to data to file, default to be 5.
"""
function GradientDescent(
	gll :: Paras, gds :: GDStatistics, datadir :: String, prevData :: Array{Float64, 2}=zeros(Float64, 1, 0), saveRate :: Int=5
	)

	# datadir = joinpath(@__DIR__, "data", string(gll.N1)*"_"*string(gll.N2), "D_"*string(round(gll.D_ion, digits=2)), "H_"*string(round(gll.h_ext, digits=3)))
	isdir(datadir) || mkpath(datadir)

	enable_MPI = false
	commSize   = 1
	rank       = 0

	# MPI initialization, not used with a single thread on laptop
	if MPI.Initialized()
		commSize = MPI.Comm_size(MPI.COMM_WORLD)
		rank     = MPI.Comm_rank(MPI.COMM_WORLD)
		if commSize > 1
			rank == 0 && @printf("MPI detected. Random + global minimization across %d simulations.\n", commSize)
			enable_MPI = true
		end
	end

	# init the opt energy and angles
	opt_energy = 1e6
	opt_angles = zeros(Float64, 4*gll.Ntot)
	# NLopt set up
	opt = Opt(:LD_LBFGS, 4*gll.Ntot)
	opt.ftol_rel = 0.0
	opt.ftol_abs = 0.0
	opt.min_objective = (xx, GG) -> Energy_obj(xx, GG, gll)
	# init the init angle array
	init_angles = zeros(Float64, 4*gll.Ntot)

	# minimization using previous data, if prevData is empty, the calculation is
	# not performed
	lenPrevData = length(prevData[1, :])
	if lenPrevData >= 1
		println("Optimization using previous data on rank ", rank, ".")
		Pstep = 1
		while Pstep <= lenPrevData
			println("Prev run ", Pstep, "-th in rank ", rank)
			init_angles = prevData[:, Pstep]
			(minf, minx, ret) = optimize(opt, init_angles)
			if minf < opt_energy
				println("Find better solution using previous results in step ", Pstep, " with energy per site= ", minf/gll.Ntot)
				opt_energy = minf
				opt_angles = minx
				ret_global = ret
				gds.numPacc += 1
			end
			Pstep += 1
		end
	end

	Iter = 1

	while Iter <= gds.numIter
		# minimization with random initial conditions within the same rank
		# this is used on a laptop.
		println("Iteration... ", Iter)
		Rstep = 1
		println("Optimization with random initial conditions in rank ", rank)
		while Rstep <= gds.numRand
			println("Rand run ", Rstep, "-th in rank ", rank)
			for site = 1:gll.Ntot
				# uniform sampling on CP^{2}
				rs = rand(Float64, 4)
				init_angles[(site-1)*4+1] = asin(rs[1]^(0.25))
				init_angles[(site-1)*4+2] = asin(rs[2]^(0.5))
				init_angles[(site-1)*4+3] = 2.0 * pi * rs[3]
				init_angles[(site-1)*4+4] = 2.0 * pi * rs[4]
			end
			(minf, minx, ret) = optimize(opt, init_angles)
			if minf < opt_energy
				println("Find new local minimum in random test ", Rstep, " with energy per site= ", minf/gll.Ntot)
				opt_energy = minf
				opt_angles = minx
				ret_global = ret
				gds.numRacc += 1
			end
			Rstep % saveRate == 0 && WriteResults(gll, gds, opt_angles, datadir, Iter)
			Rstep += 1
		end
		if enable_MPI
			opt_angles_array = zeros(Float64, 4*gll.Ntot, commSize)
			opt_angles_array[:, rank+1] = opt_angles
			MPI.Allgather!(UBuffer(opt_angles_array, 4*gll.Ntot), MPI.COMM_WORLD)
			Gstep = 1
			println("Optimization with `global` initial conditions in rank ", rank)
			while Gstep <= commSize
				println("Global run ", Gstep, "-th in rank ", rank)
				init_angles = opt_angles_array[:, Gstep]
				(minf, minx, ret) = optimize(opt, init_angles)
				if minf < opt_energy
					println("Find new local minimum in ", Gstep, " with energy per site= ", minf/gll.Ntot)
					opt_energy = minf
					opt_angles = minx
					ret_global = ret
					gds.numGacc += 1
				end
				Gstep % saveRate == 0 && WriteResults(gll, gds, opt_angles, datadir, Iter)
				Gstep += 1
			end
		end

		# O_vals = zeros(Float64, 4*gll.Ntot)
		# O_site = zeros(Float64, 4)

		# for ii = 1:gll.Ntot
			# θ  = opt_angles[(ii-1)*4+1]
			# φ  = opt_angles[(ii-1)*4+2]
			# α₁ = opt_angles[(ii-1)*4+3]
			# α₂ = opt_angles[(ii-1)*4+4]
			# O_ave4!(θ, φ, α₁, α₂, O_site)
			# O_vals[(ii-1)*4+1:(ii-1)*4+4] = O_site
		# end

		WriteResults(gll, gds, opt_angles, datadir, Iter)

		Iter += 1
	end

end
