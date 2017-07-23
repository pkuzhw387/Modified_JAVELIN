import numpy as np 

yr = 365.25
tobs = 10.0

def covmatpmap_bit(mat, jd1, jd2, id1, id2, sigma, tau, slagarr, swidarr, scalearr, nx, ny, cmin, cmax, symm, model='DRW'):

	if cmax == -1:
		cmax = ny

	scale_hidden = scalearr[2]
	# for line ids(id>=2), the corresponding scale_hidden value is scale_hidden[id - 2]
	# scale_hidden = [scalearr[2 * k] for k in range(1, (len(scalearr) - 1) / 2 + 1)]
	# print "scale_hidden is ", scale_hidden

	if symm:
		if model == 'DRW':
			for j in range(cmin, cmax):
				slag2 = slagarr[id2[j] - 1]
				swid2 = swidarr[id2[j] - 1]
				scale2 = scalearr[id2[j] - 1]

				# for i in range(0, nx):
				for i in range(0, j + 1):
					slag1 = slagarr[id1[i] - 1]
					swid1 = swidarr[id1[i] - 1]
					scale1 = scalearr[id1[i] - 1]
					mat[i, j] = covmatpmapij_DRW(mat[i, j], id1[i], id2[j], jd1[i], jd2[j], sigma, tau, \
							 		slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden)
					mat[j, i] = mat[i, j]


		elif model == 'pow-law':
			for j in range(cmin, cmax):
				slag2 = slagarr[id2[j] - 1]
				swid2 = swidarr[id2[j] - 1]
				scale2 = scalearr[id2[j] - 1]

				for i in range(0, nx):
				# for i in range(0, j + 1):
					slag1 = slagarr[id1[i] - 1]
					swid1 = swidarr[id1[i] - 1]
					scale1 = scalearr[id1[i] - 1]
					covmatpmapij_PL(mat[i, j], id1[i], id2[j], jd1[i], jd2[j], A, gamma, \
					 				slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden)


	else:
		for j in range(cmin, cmax):
			slag2 = slagarr[id2[j] - 1]
			swid2 = swidarr[id2[j] - 1]
			scale2 = scalearr[id2[j] - 1]

			for i in range(0, nx):
				slag1 = slagarr[id1[i] - 1]
				swid1 = swidarr[id1[i] - 1]
				scale1 = scalearr[id1[i] - 1]
				if model == 'DRW':
					covmatpmapij_DRW(mat[i, j], id1[i], id2[j], jd1[i], jd2[j], sigma, tau, \
							 	slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden)
				elif model == 'pow-law':
					A = sigma
					gamma = tau
					covmatpmapij_PL(mat[i, j], id1[i], id2[j], jd1[i], jd2[j], A, gamma, \
							 	slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden)

def covmatpmapij_DRW(covij, id1, id2, jd1, jd2, sigma, tau, slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden):
	imax = np.max([id1, id2])
	imin = np.min([id1, id2])
	if imin <= 0.0:
		print "ids cannot be smaller than 1!"
		covij = -1.0
		return

	if imin == imax:
		if imin == 1:
			covij = np.exp(-abs(jd1 - jd2) / tau)
		else:
			covij = scale_hidden * scale_hidden * np.exp(-abs(jd1 - jd2) / tau) +\
        			scale_hidden * scale2 * np.exp(-abs(jd1 - jd2 + slag2) / tau) +\
        			scale_hidden * scale1 * np.exp(-abs(jd1 - jd2 - slag1) / tau) +\
        			scale1 * scale2 * np.exp(-abs(jd1 - jd2) / tau)

	else:
		covij = scale_hidden * np.exp(-abs(jd1 - jd2) / tau)
		if id1 == 1 and id2 >= 2:
			covij = covij + scale2 * np.exp(-abs(jd1 - jd2 + slag2) / tau)
		elif id2 == 1 and id1 >= 2:
			covij = covij + scale1 * np.exp(-abs(jd2 - jd1 + slag1) / tau)
	covij = 0.5 * sigma**2 * covij
	return covij

# # exact the same version as the fortran code
# def covmatpmapij_DRW(covij, id1, id2, jd1, jd2, sigma, tau, slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden):
# 	imax = np.max([id1, id2])
# 	imin = np.min([id1, id2])
# 	if imin <= 0.0:
# 		print "ids cannot be smaller than 1!"
# 		covij = -1.0
# 		return

# 	if imin == imax:
# 		if imin == 1:
# 			covij = getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, slag2, scale2)
# 		else:
# 			covij = getcmat_delta_DRW(id1, id2, jd1, jd2, tau, 0.0, scale_hidden, 0.0, scale_hidden)
# 			if swid2 <= 0.01:
# 				covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, 0.0, scale_hidden, slag2, scale2)
# 			else:
# 				covij += getcmat_lc_DRW(id1, id2, jd1, jd2, tau, 0.0, 0.0, scale_hidden, slag2, swid2, scale2)

# 			if swid1 <= 0.01:
# 				covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, 0.0, scale_hidden)
# 			else:
# 				covij += getcmat_lc_DRW(id1, id2, jd1, jd2, tau, slag1, swid1, scale1, 0.0, 0.0, scale_hidden)
	
# 			if swid1 <= 0.01:
# 				covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, slag2, scale2)
# 			else:
# 				covij += getcmat_lauto_DRW(id1, jd1, jd2, tau, slag1, swid1, scale1)
# 	else:
# 		covij = getcmat_delta_DRW(id1, id2, jd1, jd2, tau, 0.0, 1.0, 0.0, scale_hidden)
# 		twidth = np.max([swid1, swid2])
# 		if twidth <= 0.01:
# 			covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, slag2, scale2)
# 		else:
# 			covij += getcmat_lc_DRW(id1, id2, jd1, jd2, tau, slag1, swid1, scale1, slag2, swid2, scale2)
# 	covij *= 0.5 * sigma**2
# 	return covij

# def covmatpmapij_DRW(covij, id1, id2, jd1, jd2, sigma, tau, slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden):
# 	imax = np.max([id1, id2])
# 	imin = np.min([id1, id2])

# 	if imin <= 0:
# 		print "ids cannot be smaller than 1!"
# 		covij = -1.0
# 		return covij

# 	if imin == 1:	
# 		if imax == 1:
# 			# id1 = id2 = 1
# 			# i.e. continuum band auto correlation: cov(c0i, c0j)
# 			# pure delta function version
# 			# covij = np.exp(-abs(jd1 - jd2) / tau)

# 			# realistic version
# 			covij = getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, slag2, scale2)


# 		else:
# 			# between two epochs of continuum and one of the lines
# 			# covij = scale_hidden * np.exp(-abs(jd1 - jd2) / tau)
# 			# In multi-line version, the scale_hidden is a array, so which scale_hidden to use
# 			#  must be specified  and the index rule is given in the covmatpmap function.
# 			twidth = np.max([swid1, swid2])
# 			if id1 == 1 and id2 >= 2:
# 				# pure delta version
# 				# covij = scale_hidden[id2-2] * np.exp(-abs(jd1 - jd2) / tau) + scale2 * np.exp(-abs(jd1 - jd2 + slag2) / tau)
				
# 				# realistic version
# 				covij = getcmat_delta_DRW(id1, id2, jd1, jd2, tau, 0.0, 1.0, 0.0, scale_hidden[id2-2])

# 			elif id2 == 1 and id1 >= 2:
# 				# pure delta version
# 				# covij = scale_hidden[id1-2] * np.exp(-abs(jd1 - jd2) / tau) + scale1 * np.exp(-abs(jd2 - jd1 + slag1) / tau)

# 				# realistic version
# 				covij = getcmat_delta_DRW(id1, id2, jd1, jd2, tau, 0.0, scale_hidden[id1-2], 0.0, 1.0)

# 			if twidth <= 0.01:
# 				covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, slag2, scale2)
# 			else:
# 				covij += getcmat_lc_DRW(id1, id2, jd1, jd2, tau, slag1, swid1, scale1, slag2, swid2, scale2)
# 	else:		
# 			# id1, id2 >= 2
# 			# i.e. line band auto correlation: cov(cmi,cnj) + cov(cmi,lnj) + cov(lmi,cnj) + cov(lmi,lnj)
			
# 			# pure delta function version
# 			# covij = scale_hidden[id1-2]**2 * np.exp(-abs(jd1 - jd2) / tau) + \
# 			# 		scale_hidden[id1-2] * scale2 * np.exp(-abs(jd1 - jd2 + slag2) / tau) + \
# 			# 		scale_hidden[id1-2] * scale1 * np.exp(-abs(jd1 - jd2 - slag1) / tau) + \
# 			# 		scale1 * scale2 * np.exp(-abs(jd1 - jd2) / tau)

# 			# realistic version
# 			# cov(cmi,cnj)
# 			covij = getcmat_delta_DRW(id1, id2, jd1, jd2, tau, 0.0, scale_hidden[id1-2], 0.0, scale_hidden[id2-2])
# 			# cov(cmi,lnj)
# 			if swid2 <= 0.01:
# 				covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, 0.0, scale_hidden[id1-2], slag2, scale2)
# 			else:
# 				covij += getcmat_lc_DRW(id1, id2, jd1, jd2, tau, 0.0, 0.0, scale_hidden[id1-2], slag2, swid2, scale2)
# 			# cov(lmi,cnj)
# 			if swid1 <= 0.01:
# 				covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, 0.0, scale_hidden[id2-2])
# 			else:
# 				covij += getcmat_lc_DRW(id1, id2, jd1, jd2, tau, slag1, swid1, scale1, 0.0, 0.0, scale_hidden[id2-2])
# 			# cov(lmi,lnj)
# 			if id1 == id2:
# 				if swid1 <= 0.01:
# 					covij += getcmat_delta_DRW(id1, id2, jd1, jd2, tau, slag1, scale1, slag2, scale2)
# 				else:
# 					covij += getcmat_lauto_DRW(id1, jd1, jd2, tau, slag1, swid1, scale1)
# 			else:
# 				covij += getcmat_lcross_DRW(id1, id2, jd1, jd2, tau, slag1, swid1, scale1, slag2, swid2, scale2)
# 	covij = 0.5 * sigma**2 * covij
# 	return covij

def getcmat_delta_DRW(id1, id2, jd1, jd2, tau, tspike1, scale1, tspike2, scale2):
	getcmat_delta = np.exp(-abs(jd1 - jd2 - tspike1 + tspike2) / tau)
	getcmat_delta = abs(scale1 * scale2) * getcmat_delta
	return getcmat_delta

def getcmat_lc_DRW(id1, id2, jd1, jd2, tau, slag1, swid1, scale1, slag2, swid2, scale2):
	if id1 == 1 and id2 >= 2:
		tlow = jd2 - jd1 - slag2 - 0.5 * swid2
		thig = jd2 - jd1 - slag2 + 0.5 * swid2
		emscale = abs(scale2)
		twidth = swid2

	elif id2 == 1 and id1 >= 2:
		tlow = jd1 - jd2 - slag1 - 0.5 * swid1
		thig = jd1 - jd2 - slag1 + 0.5 * swid1
		emscale = abs(scale1)
		twidth = swid1

	elif id1 >= 2 and id2 >= 2:
		if swid1 <= 0.01:
			tlow = jd2 - (jd1 - slag1) - slag2 - 0.5 * swid2
			thig = jd2 - (jd1 - slag1) - slag2 + 0.5 * swid2
			emscale = abs(scale2 * scale1)
			twidth = swid2
		elif swid2 <= 0.01:
			tlow = jd1 - (jd2 - slag2) - slag1 - 0.5 * swid1
			thig = jd1 - (jd2 - slag2) - slag1 + 0.5 * swid1
			emscale = abs(scale2 * scale1)
			twidth = swid1

	if thig <= 0:
		getcmat_lc = np.exp(thig / tau) - np.exp(tlow / tau)
	elif tlow >= 0.0:
		getcmat_lc = np.exp(-tlow / tau) - np.exp(-thig / tau)
	else: 
		getcmat_lc = 2.0 - np.exp(tlow / tau) - np.exp(-thig / tau)

	getcmat_lc = tau * emscale / twidth * getcmat_lc
	return getcmat_lc

def getcmat_lauto_DRW(id, jd1, jd2, tau, slag, swid, scale):
	twidth = swid
	emscale = scale 

	tlow = jd1 - jd2 - twidth
	tmid = jd1 - jd2
	thig = jd1 - jd2 + twidth

	if thig <= 0.0 or tlow >= 0.0:
		getcmat_lauto = np.exp(-abs(tmid) / tau) * (np.exp(0.5 * twidth / tau) - np.exp(-0.5 * twidth / tau))**2
	else:
		getcmat_lauto = -2.0 * np.exp(-abs(tmid) / tau) + np.exp(-twidth / tau) * (np.exp(-tmid / tau) + np.exp(tmid / tau))
		if tmid >= 0.0:
			getcmat_lauto = getcmat_lauto + 2.0 * (twidth - tmid) / tau
		else:
			getcmat_lauto = getcmat_lauto + 2.0 * (twidth + tmid) / tau
	getcmat_lauto = ((emscale * tau / twidth)**2) * getcmat_lauto
	return getcmat_lauto

def getcmat_lcross_DRW(id1, id2, jd1, jd2, tau, slag1, swid1, scale1, slag2, swid2, scale2):
	twidth1 = swid1
	twidth2 = swid2

	if twidth1 >= twidth2:
		t1 = slag1 - 0.5 * twidth1
		t2 = slag1 + 0.5 * twidth1
		t3 = slag2 - 0.5 * twidth2
		t4 = slag2 + 0.5 * twidth2
		bottleneck = twidth2
		ti = jd1
		tj = jd2 
	else:
		t1 = slag2 - 0.5 * twidth2
		t2 = slag2 + 0.5 * twidth2
		t3 = slag1 - 0.5 * twidth1
		t4 = slag1 + 0.5 * twidth1
		bottleneck = twidth1
		ti = jd2 
		tj = jd1

	tlow = (ti - tj) - (t2 - t3)
	tmid1 = (ti - tj) - (t2 - t4)
	tmid2 = (ti - tj) - (t1 - t3)
	thig = (ti - tj) - (t1 - t4)

	if thig <= 0.0 or tlow >= 0.0:
		getcmat_lcross = np.exp(-abs(tlow) / tau) + np.exp(-abs(thig) / tau) \
						- np.exp(-abs(tmid1) / tau) - np.exp(-abs(tmid2)/ tau)
	else:
		getcmat_lcross = np.exp(tlow / tau) + np.exp(-thig / tau) \
						- np.exp(-abs(tmid1) / tau) - np.exp(-abs(tmid2) / tau)

		if tmid2 <= 0.0:
			getcmat_lcross = getcmat_lcross + 2 * thig / tau
		elif tmid1 <= 0.0:
			getcmat_lcross = getcmat_lcross + 2 * bottleneck / tau
		elif tlow < 0.0:
			getcmat_lcross = getcmat_lcross - 2 * tlow / tau

	getcmat_lcross = tau * tau * scale1 * scale2 / (twidth1 * twidth2) * getcmat_lcross
	return getcmat_lcross



def covmatpmapij_PL(covij, id1, id2, jd1, jd2, A, gamma, slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden):
	imax = np.max([id1, id2])
	imin = np.min([id1, id2])

	if imin <= 0:
		print "ids cannot be smaller than 1!"
		covij = -1.0
		return covij

	if imin == imax:
		if imin == 1:
			covij = tobs**gamma - 0.5 * (abs(jd1 - jd2) / yr)**gamma
		else:
			covij = scale_hidden * scale_hidden * (tobs**gamma - 0.5 * (abs(jd1 - jd2) / yr)**gamma) +\
         			scale_hidden * scale2 * (tobs**gamma - 0.5 * (abs(jd1 - jd2 + slag2) / yr)**gamma) +\
         			scale_hidden * scale1 * (tobs**gamma - 0.5 * (abs(jd1 - jd2 - slag1) / yr)**gamma) +\
         			scale1 * scale2 *(tobs**gamma - 0.5 * (abs(jd1 - jd2) / yr)**gamma)

	else:
		covij = scale_hidden * (tobs**gamma - 0.5 * (abs(jd1 - jd2) / yr)**gamma)
		if id1 == 1 and id2 >= 2:
			covij = covij + scale2 * (tobs**gamma - 0.5 * (abs(jd1 - jd2 + slag2) / yr)**gamma)
		elif id2 == 1 and id1 >= 2:
			covij = covij + scale1 * (tobs**gamma - 0.5 * (abs(jd2 - jd1 + slag1) / yr)**gamma)
		covij = A**2 * covij
		return covij




# def covmatpmapij_PL(covij, id1, id2, jd1, jd2, A, gamma, slag1, swid1, scale1, slag2, swid2, scale2, scale_hidden):
# 	imax = np.max([id1, id2])
# 	imin = np.min([id1, id2])

# 	if imin <= 0:
# 		print "ids cannot be smaller than 1!"
# 		covij = -1.0
# 		return covij

# 	if imin == 1:	
# 		if imax == 1:
# 			# id1 = id2 = 1
# 			# i.e. continuum band auto correlation: cov(c0i, c0j)
# 			# pure delta function version
# 			# covij = np.exp(-abs(jd1 - jd2) / tau)

# 			# realistic version
# 			covij = getcmat_delta_PL(id1, id2, jd1, jd2, gamma, slag1, scale1, slag2, scale2)


# 		else:
# 			# print "id1 and id2: ", id1, id2
# 			# between two epochs of continuum and one of the lines
# 			# covij = scale_hidden * np.exp(-abs(jd1 - jd2) / tau)
# 			# In multi-line version, the scale_hidden is a array, so which scale_hidden to use
# 			#  must be specified  and the index rule is given in the covmatpmap function.
# 			twidth = np.max([swid1, swid2])
# 			if id1 == 1 and id2 >= 2:
# 				# pure delta version
# 				# covij = scale_hidden[id2-2] * np.exp(-abs(jd1 - jd2) / tau) + scale2 * np.exp(-abs(jd1 - jd2 + slag2) / tau)
				
# 				# realistic version
# 				covij = getcmat_delta_PL(id1, id2, jd1, jd2, gamma, 0.0, 1.0, 0.0, scale_hidden[id2-2])
# 				# print "covij: ", covij
# 			elif id2 == 1 and id1 >= 2:
# 				# pure delta version
# 				# covij = scale_hidden[id1-2] * np.exp(-abs(jd1 - jd2) / tau) + scale1 * np.exp(-abs(jd2 - jd1 + slag1) / tau)

# 				# realistic version
# 				covij = getcmat_delta_PL(id1, id2, jd1, jd2, gamma, 0.0, scale_hidden[id1-2], 0.0, 1.0)

# 			if twidth <= 0.01:
# 				covij += getcmat_delta_PL(id1, id2, jd1, jd2, gamma, slag1, scale1, slag2, scale2)
# 			else:
# 				# print covij
# 				covij += getcmat_lc_PL(id1, id2, jd1, jd2, gamma, slag1, swid1, scale1, slag2, swid2, scale2)
# 	else:		
# 			# id1, id2 >= 2
# 			# i.e. line band auto correlation: cov(cmi,cnj) + cov(cmi,lnj) + cov(lmi,cnj) + cov(lmi,lnj)
			
# 			# pure delta function version
# 			# covij = scale_hidden[id1-2]**2 * np.exp(-abs(jd1 - jd2) / tau) + \
# 			# 		scale_hidden[id1-2] * scale2 * np.exp(-abs(jd1 - jd2 + slag2) / tau) + \
# 			# 		scale_hidden[id1-2] * scale1 * np.exp(-abs(jd1 - jd2 - slag1) / tau) + \
# 			# 		scale1 * scale2 * np.exp(-abs(jd1 - jd2) / tau)

# 			# realistic version
# 			# cov(cmi,cnj)
# 			covij = getcmat_delta_PL(id1, id2, jd1, jd2, gamma, 0.0, scale_hidden[id1-2], 0.0, scale_hidden[id2-2])
# 			# cov(cmi,lnj)
# 			if swid2 <= 0.01:
# 				covij += getcmat_delta_PL(id1, id2, jd1, jd2, gamma, 0.0, scale_hidden[id1-2], slag2, scale2)
# 			else:
# 				covij += getcmat_lc_PL(id1, id2, jd1, jd2, gamma, 0.0, 0.0, scale_hidden[id1-2], slag2, swid2, scale2)
# 			# cov(lmi,cnj)
# 			if swid1 <= 0.01:
# 				covij += getcmat_delta_PL(id1, id2, jd1, jd2, gamma, slag1, scale1, 0.0, scale_hidden[id2-2])
# 			else:
# 				covij += getcmat_lc_PL(id1, id2, jd1, jd2, gamma, slag1, swid1, scale1, 0.0, 0.0, scale_hidden[id2-2])
# 			# cov(lmi,lnj)
# 			if id1 == id2:
# 				if swid1 <= 0.01:
# 					covij += getcmat_delta_PL(id1, id2, jd1, jd2, gamma, slag1, scale1, slag2, scale2)
# 				else:
# 					covij += getcmat_lauto_PL(id1, jd1, jd2, gamma, slag1, swid1, scale1)
# 			else:
# 				covij += getcmat_lcross_PL(id1, id2, jd1, jd2, gamma, slag1, swid1, scale1, slag2, swid2, scale2)
# 	covij = A**2 * covij
# 	return covij

def getcmat_delta_PL(id1, id2, jd1, jd2, gamma, tspike1, scale1, tspike2, scale2):
	getcmat_delta = tobs**gamma - 0.5 * abs((jd1 - jd2 - tspike1 + tspike2) / yr)**gamma
	getcmat_delta = abs(scale1 * scale2) * getcmat_delta
	return getcmat_delta

def getcmat_lc_PL(id1, id2, jd1, jd2, gamma, slag1, swid1, scale1, slag2, swid2, scale2):
	if id1 == 1 and id2 >= 2:
		tlow = jd2 - jd1 - slag2 - 0.5 * swid2
		thig = jd2 - jd1 - slag2 + 0.5 * swid2
		emscale = abs(scale2)
		twidth = swid2

	elif id2 == 1 and id1 >= 2:
		tlow = jd1 - jd2 - slag1 - 0.5 * swid1
		thig = jd1 - jd2 - slag1 + 0.5 * swid1
		emscale = abs(scale1)
		twidth = swid1

	elif id1 >= 2 and id2 >= 2:
		if swid1 <= 0.01:
			tlow = jd2 - (jd1 - slag1) - slag2 - 0.5 * swid2
			thig = jd2 - (jd1 - slag1) - slag2 + 0.5 * swid2
			emscale = abs(scale2 * scale1)
			twidth = swid2

		elif swid2 <= 0.01:
			tlow = jd1 - (jd2 - slag2) - slag1 - 0.5 * swid1
			thig = jd1 - (jd2 - slag2) - slag1 + 0.5 * swid1
			emscale = abs(scale2 * scale1)
			twidth = swid1

	if thig <= 0:
		getcmat_lc = (-thig)**(gamma + 1) - (-tlow)**(gamma + 1) + 2 * thig * tobs**gamma \
					+ 2 * gamma * thig * tobs**gamma - 2 * tlow * tobs**gamma - 2 * gamma * tlow * tobs**gamma
	elif tlow >= 0.0:
		getcmat_lc = -thig**(gamma + 1) + tlow**(gamma + 1) + 2 * thig * tobs**gamma \
					+ 2 * gamma * thig * tobs**gamma - 2 * tlow * tobs**gamma - 2 * gamma * tlow * tobs**gamma
	else: 
		getcmat_lc = -thig ** (1 + gamma) - (-tlow)**(1 + gamma) + 2 * thig * tobs**gamma \
					+ 2 * gamma * thig * tobs**gamma - 2 * tlow * tobs**gamma - 2 * gamma * tlow * tobs**gamma
	getcmat_lc = emscale / twidth * getcmat_lc / (2 + 2 * gamma)
	return getcmat_lc

def getcmat_lauto_PL(id, jd1, jd2, gamma, slag, swid, scale):
	twidth = swid
	emscale = scale 
	# print "gamma is: ", gamma
	tlow = jd1 - jd2 - twidth
	tmid = jd1 - jd2 
	thig = jd1 - jd2 + twidth
	if thig <= 0.0:
		getcmat_lauto = (-tlow)**(gamma + 2) + tlow**2 * (tlow - 2 * tmid)**gamma \
						- 2 * tlow * (2 * (tlow - 2 * tmid)**gamma * tmid + (2 + 3 * gamma + gamma**2) * (-twidth) * tobs**gamma)\
						+ 2 * tmid * ((-tmid)**(gamma + 1) + 2 * (tlow - 2 * tmid)**gamma * tmid + (2 + 3 * gamma + gamma**2) * (-twidth) * tobs**gamma)
		getcmat_lauto *= -1.0
	elif tlow >= 0:
		getcmat_lauto = -tlow ** (gamma + 2) + 4 * tlow * tmid * ((-tlow + 2 * tmid)**gamma - (2 + 3 * gamma + gamma**2) * tobs**gamma) \
						+ 2 * tmid**2 * (tmid**gamma - 2 * (-tlow + 2 * tmid)**gamma + (2 + 3 * gamma + gamma**2) * tobs**gamma) \
						+ tlow**2 (-(-tlow + 2 * tmid)**gamma + 2 * (2 + 3 * gamma + gamma**2) * tobs**gamma)
	else:
		if tmid <= 0.0:
			getcmat_lauto = -tlow**2 * ((-tlow)**gamma + (-tlow + 2 * tmid)**gamma - 2 * (2 + 3 * gamma + gamma**2) * tobs**gamma)\
							+ 4 * tlow * tmid * ((-tlow + 2 * tmid)**gamma - (2 + 3 * gamma + gamma**2) * tobs**gamma)\
							+ 2 * tmid**2 * ((-tmid)**gamma - 2 * (-tlow + 2 * tmid)**gamma + (2 + 3 * gamma + gamma**2) * tobs**gamma)
		else:
			getcmat_lauto = -tlow**2((-tlow)**gamma + (-tlow + 2 * tmid)**gamma - 2 * (2 + 3 * gamma + gamma**2) * tobs**gamma)\
							+ 4 * tlow * tmid * ((-tlow + 2 * tmid)**gamma - (2 + 3 * gamma + gamma**2) * tobs**gamma)\
							+ 2 * tmid**2 * ((tmid)**gamma - 2 * (-tlow + 2 * tmid)**gamma + (2 + 3 * gamma + gamma**2) * tobs**gamma)

	getcmat_lauto = ((emscale / twidth)**2) * getcmat_lauto / ((2 * gamma + 2) * (gamma + 2))
	return getcmat_lauto


def getcmat_lcross_PL(id1, id2, jd1, jd2, gamma, slag1, swid1, scale1, slag2, swid2, scale2):
	twidth1 = swid1
	twidth2 = swid2

	if twidth1 >= twidth2:
		t1 = slag1 - 0.5 * twidth1
		t2 = slag1 + 0.5 * twidth1
		t3 = slag2 - 0.5 * twidth2
		t4 = slag2 + 0.5 * twidth2
		bottleneck = twidth2
		ti = jd1
		tj = jd2 
	else:
		t1 = slag2 - 0.5 * twidth2
		t2 = slag2 + 0.5 * twidth2
		t3 = slag1 - 0.5 * twidth1
		t4 = slag1 + 0.5 * twidth1
		bottleneck = twidth1
		ti = jd2 
		tj = jd1

	tlow = (ti - tj) - (t2 - t3)
	tmid1 = (ti - tj) - (t2 - t4)
	tmid2 = (ti - tj) - (t1 - t3)
	thig = (ti - tj) - (t1 - t4)