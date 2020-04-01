def calc_Qrad(Aio, Aoi):
	C = 5.67/(1/epsio - 1 + 1/phi + (1/epsoi-1)*Aio/Aoi)
	return C
