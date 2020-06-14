#To find
#1. Subroutine to calculate Twio and Twoi, Twii, Tm
#2. How to calculate Pin and Psv
#3. How to find velocity (c and cout)
#4. How to find rho and mu
#5. 2nd and 3rd term in equation 36
#6. deltaT in equation 22

from pressure_drop import *
from Heat_transfer import *

def calc_area(d):
	return np.pi*(d**2)/4
def calc_larea(d):
	return np.pi*d*l
def calc_phi2(x,mul,mug,rhol,rhog):
	return (1 + x*(mul - mug)/mug)**(-0.25) * (1 + x*(rhol/rhog - 1))
def dif(TMLI, Aio, Aoi, rio, roi, Twio, spacer, dio, Twoi, NuLO, Twii, dii, Tb):
	return calc_Qrad(Aio, Aoi, rio, roi, Twio, TMLI) + calc_Qspacer(spacer, dio, Twio, Twoi) - calc_QconvA(NuLO, Twii, dii, Tb)


#Initial values
psv =
hLI = 
mHetot = 
ltot = 
TLHein =
Pin =
cin =
rhoin = 
dtubeout =
dtubein = 
Atubeout = calc_area(dtubeout)
Atubein = calc_area(dtubein)

n=0
Pinj = np.zeros(integer(ltot/l))
TLHeinj = np.zeros(integer(ltot/l))
hinj = np.zeros(integer(ltot/l))
lz = np.zeros(integer(ltot/l))
TLHeinj[0] = TLHein
hinj[0] = hLI
Pinj[0] = Pin

while True:
	j = 0
	mtot = rhoin * cin * Atubein	
	xin = 1
	c = cin
	Pfr = 0
	Pfit = 0
	while (j<integer(ltot/l)):
		K = 0 #no consideration of fitting loss
		spacer = 0 #no consideration of spacer loss
		#calculate K and spacer by the following if conditions if want to consider		
		if lz< :
			tube = 's'
			dii =
			gamma = 
			Twio = 4 #Assuming no gradient
			dio = 
			doi = 
			AMLIo = 
		elif lz<:
			tube = 'r'
			dii =
			gamma =
			dio = 
			doi = 
			Twio = 300
			AMLIo = 
		else:
			tube = 's'
			dii =
			Twio = 4
			gamma =
			dio = 
			doi = 
			AMLIo = 
		Aii = calc_area(dii)
		Alii= calc_larea(dii)
		Aio = calc_area(dio)			
		Aoi = calc_area(doi)
		rhol =
		mul = 
		Tb = 
		if (TLHeinj[j] < Tb):
			x = 0
			phisq = 1
			rhog = NA
			mug = NA
		else:
			x = xin
			mug = 
			rhog =
			phisq = calc_phi2(x,mul,mug,rhol,rhog)
		Poutk = calc_Pout(Pinj[j], mtot, tube, dii, x, xin, c, K, rhol,rhog,phisq,gamma,mul)) #Pacc will be zero for this
		Re = calc_rhotp(x,rhol,rhog)*c*dii/calc_mutp(x,mul,mug)
		cp =
		cpl =
		k =
		lamda = 
		Pr = calc_mutp(x,mul,mug)*cp/k
		Prl = mul*cpl/k
		hsp = calc_hsp(Re,Pr,lamda,dii)
		NuLO = hsp*dii/k
		hl = 
		hg =
		r = hg - hl	
		while True: 	
			Poutk1 = Poutk		
			TMLI = Tb + 0.5
			while True:
				TMLI = optimize.root_scalar(dif, args = (Aio, Aoi, dio/2, doi/2, Twio, spacer, dio, Twoi, NuLO, Twii, dii, Tb))				
				if calc_QONBA(deltaT,NuLO,r,dii,x) < calc_QSNBA(NuLO, Twii, dii, cpl,Tb,Prl,r):
					htp = calc_hFCE(x, NuLO, dii)
				else:
					htp = calc_hSNB(NuLO, Twii, dii, cpl,Tb,Prl,r)
				alphasp = calc_alpha(1,rhol,rhog)
				alphatp = calc_alpha(x,rhol,rhog)	
				if alphasp>alphatp:
					alpha =	alphasp
				else:
					alpha =	alphatp		
				RMLI = (TMLI-Twoi)/(calc_QMLIrad(AMLIo, Tm, TMLI, Twoi) + calc_QMLIcond(doi, TMLI, Twoi))
				Rtube = calc_Rtube(dii/2, doi/2, htp)
				Qrad = calc_Qrad(Aio, Aoi, dio/2, doi/2, Twio, TMLI)
				Qspacer = calc_Qspacer(spacer, dio, Twio, Twoi)
				Qconv = calc_QconvA(NuLO, Twii, dii, Tb) * Alii
				if (Qconv - Qrad - Qspacer < 10**(-12)): #checking equality
					break
			Rtube = 1/(2*np.pi*l)*(1/(calc_alpha(x,rhol,rhog)*dii/2) + 1/) + 1/lambdatube*np.log(doi/dii)
			RMLI = (TMLI-Twoi)/(calc_QMLIrad(AMLIo, Tm, TMLI, Twoi) + calc_QMLIcond(doi, TMLI, Twoi))
			Qtot = (TMLI-Tb)/(Rtube+RMLI) + calc_Qspacer(spacer, dio, Twio, Twoi) + Qvalve
			mtot = ((1-x)*rhol + x*rhog)*c*Aii
			hout = hin[j] + Qtot/mtot -0.5*(cout**2-c**2)- g*np.sin(gamma)*l
			xout = (hout - hl)/(hg - hl)
			Poutk =  calc_Pout(Pinj[j], mtot, tube, dii, xout, xin, c, K, rhol,rhog,phisq,gamma,mul))
		if (abs(Poutk - Poutk1) < 10**(-12)): #in MPa
			break
		j = j+1
		pfr = pfr + calc_Pfr(tube, dii, x, c, rhol, rhog, mul, phisq)
		print(j)
		lz[j] = lz[j-1] + l
		hin[j] = hout
		Pinj[j] = Poutk
		TLHeinj[j] = #calculate using hin and Pinj
		
		xin = x
	rhoout = calc_rhotp(x,rhol,rhog)	
	k = k+1

	coutn1 = mHetot/(rhoout*Atubeout))
	mRG = x*rhog*Atubeout*coutn1
	ptv = 0.0975*np.exp(59.59*mRG) =
	n = n+1
	print(n)
	coutn = np.sqrt((pin + rhol*g*hLI + rhog*g*hLg - pout - rhoout*g*hsv - pfr - pfit)*2/rhoout)
	mHetot = coutn*rhoout*Atubeout 
	if (abs(coutn - coutn1)<0.01):
		break
#plotting values (p,h,T)
plt.plot(lz,Pinj)
