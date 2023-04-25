function getMagnetization(O_vals::Vector{Float64}, gll::Paras)
	Ox, Oy, Oz = 0.0, 0.0, 0.0
	for ii ∈ 1:gll.Ntot
		O1, O2, O3 = O_vals[(ii-1)*4+1], O_vals[(ii-1)*4+2], O_vals[(ii-1)*4+3]
		Ox += O1
		Oy += O2
		Oz += O3
	end
	sqrt(Ox^2+Oy^2+Oz^2) / gll.Ntot
end

function getEnergy(O_vals::Vector{Float64}, gll::Paras)
	Energy = 0.0
	for ii ∈ 1:gll.Ntot
		@views begin
			site_nei = gll.Nei_table[:, ii]
			O1, O2, O3, O8 = O_vals[(ii-1)*4+1], O_vals[(ii-1)*4+2],
			O_vals[(ii-1)*4+3], O_vals[(ii-1)*4+4]
		end
		for nn = 1:gll.N_bond
			sitep = site_nei[nn]
			@views begin
				mat = gll.J_mat[:, :, nn]
				O1p, O2p, O3p = O_vals[(sitep-1)*4+1], O_vals[(sitep-1)*4+2],
				O_vals[(sitep-1)*4+3]
			end
			# insert 0.5 to avoid over-counting.
			# each bond is counted twice
			Energy += 0.5 * ([O1 O2 O3] * mat * [O1p; O2p; O3p])[1]
		end
		Energy += gll.D_ion/sqrt(3.0) * (O8 + 2.0*sqrt(3.0)/3.0) - gll.h_ext * O3
	end
	Energy / gll.Ntot
end
