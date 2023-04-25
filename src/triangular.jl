function Jmat_init(
	J_intra  :: Vector{Float64},
	max_dist :: Int,
	lat      :: Triangular
	)

	Num_intra = [6, 6, 6]
	Num_previ = [0, 6, 12]
	Num_tot   = [6, 12, 18]
	Jmat = zeros(Float64, 3, 3, Num_tot[max_dist])

	# C3
	sym = [-0.5 -0.5*sqrt(3.0) 0.0; 0.5*sqrt(3.0) -0.5 0.0; 0.0 0.0 1.0]
	tsym = transpose(sym)

	for nn = 1:max_dist
		J, Delta_y, Delta_z, J_xy, J_xz, J_yz =
		J_intra[(nn-1)*6+1], J_intra[(nn-1)*6+2],
		J_intra[(nn-1)*6+3], J_intra[(nn-1)*6+4],
		J_intra[(nn-1)*6+5], J_intra[(nn-1)*6+6]
		tmp1 = [J J_xy J_xz; J_xy J*Delta_y J_yz; J_xz J_yz J*Delta_z]
		tmp2 = sym * tmp1 * tsym
		tmp3 = sym * tmp2 * tsym

		# First in-plane bond of the `nn`-nearest neighbor is chosen be
		# parallel to the x-axis, then rotate counter-clockwise:
		# applying the C6 operator
		Jmat[:, :, Num_previ[nn]+1] = tmp1
		Jmat[:, :, Num_previ[nn]+4] = tmp1
		Jmat[:, :, Num_previ[nn]+2] = tmp3
		Jmat[:, :, Num_previ[nn]+5] = tmp3
		Jmat[:, :, Num_previ[nn]+3] = tmp2
		Jmat[:, :, Num_previ[nn]+6] = tmp2

	end

	Jmat
end

function neighbor_idx_table(
	N1 :: Int,
	N2 :: Int,
	max_dist :: Int,
	lat   :: Triangular,
	basis :: Int
	)

	# basis == 0, J1-J2 model
	# basis == 1, J1-J3 model
	@assert basis == 0 || basis == 1

	function coor2site(a, b)
		idx = a + (b-1) * N1
	end

	function site2coor(idx)
		b = fld(idx-1, N1) + 1
		a = mod(idx-1, N1) + 1
		a, b
	end

	Ntot = N1 * N2
	Num_tot   = [6, 12, 18]
	idx_table = zeros(Int, Num_tot[max_dist], Ntot)

	for site = 1:Ntot
		ii, jj = site2coor(site)
		# nn1
		iip1 = mod(ii, N1) + 1
		jjp1 = mod(jj, N2) + 1
		iim1 = mod(ii-2, N1) + 1
		jjm1 = mod(jj-2, N2) + 1
		# nn2
		iip2 = mod(ii+1, N1) + 1
		jjp2 = mod(jj+1, N2) + 1
		iim2 = mod(ii-3, N1) + 1
		jjm2 = mod(jj-3, N2) + 1

		for nn = 1:max_dist
			if nn == 1
				# nearest neighbor
				if basis == 0
					nn1 = coor2site(iip1, jj)
					nn2 = coor2site(ii, jjp1)
					nn3 = coor2site(iim1, jjp1)
					nn4 = coor2site(iim1, jj)
					nn5 = coor2site(ii,   jjm1)
					nn6 = coor2site(iip1, jjm1)
				else
					nn1 = coor2site(iip1, jj)
					nn2 = coor2site(iip1, jjp1)
					nn3 = coor2site(ii, jjp1)
					nn4 = coor2site(iim1, jj)
					nn5 = coor2site(iim1, jjm1)
					nn6 = coor2site(ii, jjm1)
				end
				idx_table[1:6, site] = [nn1, nn2, nn3, nn4, nn5, nn6]
			elseif nn == 2
				if basis == 0
					nn1 = coor2site(iip1, jjp1)
					nn2 = coor2site(iim1, jjp2)
					nn3 = coor2site(iim2, jjp1)
					nn4 = coor2site(iim1, jjm1)
					nn5 = coor2site(iip1, jjm2)
					nn6 = coor2site(iip2, jjm1)
				else
					nn1 = coor2site(iip2, jjp1)
					nn2 = coor2site(iip1, jjp2)
					nn3 = coor2site(iim1, jjp1)
					nn4 = coor2site(iim2, jjm1)
					nn5 = coor2site(iim1, jjm2)
					nn6 = coor2site(iip1, jjm1)
				end
				idx_table[7:12, site] = [nn1, nn2, nn3, nn4, nn5, nn6]
			elseif nn == 3
				if basis == 0
					nn1 = coor2site(iip2, jj)
					nn2 = coor2site(ii,   jjp2)
					nn3 = coor2site(iim2, jjp2)
					nn4 = coor2site(iim2, jj)
					nn5 = coor2site(ii,   jjm2)
					nn6 = coor2site(iip2, jjm2)
				else
					nn1 = coor2site(iip2, jj)
					nn2 = coor2site(iip2, jjp2)
					nn3 = coor2site(ii,   jjp2)
					nn4 = coor2site(iim2, jj)
					nn5 = coor2site(iim2, jjm2)
					nn6 = coor2site(ii,   jjm2)

				end
				idx_table[13:18, site] = [nn1, nn2, nn3, nn4, nn5, nn6]

			end
		end
	end

	idx_table

end
