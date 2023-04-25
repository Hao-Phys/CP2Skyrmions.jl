"""
Returns to `O{1,2,3,8}_ave` given four angles
`theta`, `phi`, `alpha1`, and `alpha2`
"""
function O_ave4!(
	theta, phi, alpha1, alpha2,
	O_vals :: Vector{Float64}
	)

	sin2t = sin(2.0 * theta)
	cos2t = cos(2.0 * theta)
	sint_sq = sin(theta)^2

	sina1 = sin(alpha1)
	sina2 = sin(alpha2)
	cosa1 = cos(alpha1)
	cosa2 = cos(alpha2)

	sinp = sin(phi)
	cosp = cos(phi)
	cos2p = cos(2.0 * phi)

	O_vals[1] = 0.5 * sqrt(2.0) * sin2t * (cosa1*cosp + cosa2*sinp)
	O_vals[2] = 0.5 * sqrt(2.0) * sin2t * (-sina1*cosp + sina2*sinp)
	O_vals[3] = sint_sq * cos2p
	O_vals[4] = -(1.0 + 3.0*cos2t)/(2.0*sqrt(3.0))

end

function O_ave4_grad!(
	theta, phi, alpha1, alpha2,
	Grad :: Vector{Float64}
	)

	Grad[1] = sqrt(2.0) * cos(2.0*theta) * (cos(alpha1)*cos(phi) + cos(alpha2)*sin(phi))
	Grad[2] = 0.5 * sqrt(2.0) * sin(2.0*theta) * (-sin(phi)*cos(alpha1) + cos(phi)*cos(alpha2))
	Grad[3] = -0.5 * sqrt(2.0) * sin(2.0*theta) * cos(phi) * sin(alpha1)
	Grad[4] = -0.5 * sqrt(2.0) * sin(2.0*theta) * sin(phi) * sin(alpha2)
	Grad[5] = sqrt(2.0) * cos(2.0*theta) * (-sin(alpha1)*cos(phi) + sin(alpha2)*sin(phi))
	Grad[6] = 0.5 * sqrt(2.0) * sin(2.0*theta) * (sin(phi)*sin(alpha1) + cos(phi)*sin(alpha2))
	Grad[7] = -0.5 * sqrt(2.0) * sin(2.0*theta) * cos(phi) * cos(alpha1)
	Grad[8] = 0.5 * sqrt(2.0) * sin(2.0*theta) * sin(phi) * cos(alpha2)
	Grad[9] = sin(2.0*theta) * cos(2.0*phi)
	Grad[10] = -2.0 * (sin(theta))^2 * sin(2.0*phi)
	Grad[11] = 0.0
	Grad[12] = 0.0
	Grad[13] = sqrt(3.0) * sin(2.0*theta)
	Grad[14] = 0.0
	Grad[15] = 0.0
	Grad[16] = 0.0
end

function Energy_obj(
	xx   :: Vector{Float64},
	grad :: Vector{Float64},
	gll  :: Paras
	)

	O_site  = zeros(Float64, 4)
	O_grad  = zeros(Float64, 16)
	O_vals      = zeros(Float64, 4*gll.Ntot)
	O_vals_grad = zeros(Float64, 16*gll.Ntot)

	for site = 1:gll.Ntot
		theta  = xx[(site-1)*4+1]
		phi    = xx[(site-1)*4+2]
		alpha1 = xx[(site-1)*4+3]
		alpha2 = xx[(site-1)*4+4]
		O_ave4!(theta, phi, alpha1, alpha2, O_site)
		O_ave4_grad!(theta, phi, alpha1, alpha2, O_grad)
		O_vals[(site-1)*4+1:(site-1)*4+4] = O_site
		O_vals_grad[(site-1)*16+1:(site-1)*16+16] = O_grad
	end

	if length(grad) > 0
		for site = 1:gll.Ntot
			site_grad = zeros(Float64, 4)
			@views begin
				site_nei = gll.Nei_table[:, site]
				O_grad = O_vals_grad[(site-1)*16+1:(site-1)*16+16]
			end

			for nn = 1:gll.N_bond
				sitep = site_nei[nn]
				@views begin
					mat = gll.J_mat[:, :, nn]
					O1p, O2p, O3p = O_vals[(sitep-1)*4+1], O_vals[(sitep-1)*4+2],
					O_vals[(sitep-1)*4+3]
				end

				for ag = 1:4
					site_grad[ag] +=
					([O_grad[ag] O_grad[4+ag] O_grad[8+ag]] * mat *
					[O1p; O2p; O3p])[1]
				end
			end

			site_grad[1] += -gll.h_ext * O_grad[9] +
			gll.D_ion/sqrt(3.0) * O_grad[13]
			site_grad[2] += -gll.h_ext *O_grad[10]
			grad[(site-1)*4+1:(site-1)*4+4] = site_grad

		end
	end

	energy_val = 0.0

	for site = 1:gll.Ntot

		@views begin
			site_nei = gll.Nei_table[:, site]
			O1, O2, O3, O8 = O_vals[(site-1)*4+1], O_vals[(site-1)*4+2],
			O_vals[(site-1)*4+3], O_vals[(site-1)*4+4]
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
			energy_val += 0.5 * ([O1 O2 O3] * mat * [O1p; O2p; O3p])[1]
		end

		energy_val += gll.D_ion/sqrt(3.0) * (O8 + 2.0*sqrt(3.0)/3.0) - gll.h_ext * O3
	end

	return energy_val
end



"""
Computes the 8 SU(3) generators on a particular site,
given four angles.
"""
function O_ave8(theta, phi, alpha1, alpha2)
	O_ave = Vector{Float64}(undef, 8)
	O1 = 0.5 * sqrt(2.0) * sin(2.0*theta) * ( cos(alpha1)*cos(phi) + cos(alpha2)*sin(phi))
	O2 = 0.5 * sqrt(2.0) * sin(2.0*theta) * (-sin(alpha1)*cos(phi) + sin(alpha2)*sin(phi))
	O3 = sin(theta)^2 * cos(2.0*phi)
	O4 = 0.5 * sqrt(2.0) * sin(2.0*theta) * (-cos(alpha1)*cos(phi) + cos(alpha2)*sin(phi))
	O5 = 0.5 * sqrt(2.0) * sin(2.0*theta) * ( sin(alpha1)*cos(phi) + sin(alpha2)*sin(phi))
	O6 =  cos(alpha1 - alpha2) * sin(theta)^2 * sin(2.0*phi)
	O7 = -sin(alpha1 - alpha2) * sin(theta)^2 * sin(2.0*phi)
	O8 = -(1.0 + 3.0*cos(2.0*theta))/(2.0*sqrt(3.0))
	O_ave = [O1; O2; O3; O4; O5; O6; O7; O8]

end
