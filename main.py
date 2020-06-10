#To find
#1. Subroutine to calculate Twio
#2. How to calculate Pin
#3. calculation of ReLO
#3. How to find velocity
#4. How to find rho and mu


from pressure_drop import *
from Heat_transfer import *

#Initial values
psv =
hLI = 
mHetot = 
ltot = 

TLHein =
dtubeout =
def calc_area(d):
	return np.pi*(d**2)/4
def calc_phi2(x,mul,mug,rhol,rhog):
	return (1 + x*(mul - mug)/mug)**(-0.25) * (1 + x*(rhol/rhog - 1))

Pin =

n=0
Pinj = np.zeros(integer(ltot/l))
lz = np.zeros(integer(ltot/l))
Pinj[0] = Pin
hinj[0] = hLI
while True:
	j = 0
	while (j<integer(ltot/l)):
		K = 0 #no consideration of fitting loss
		spacer = 0 #no consideration of spacer loss
		#calculate K and spacer by the following if conditions if want to consider		
		if lz< :
			tube = s
			dii =
			gamma = 
		elif lz<:
			tube = r
			dii =
			gamma =
		else:
			tube = s
			dii =
			gamma =
		rhol =
		mul = 
		Tb = 
		if (TLHein < Tb):
			x = 0
			phisq = 1
			rho = rhol
			mu = mul
		else:
			x = xin
			mug = 
			rhog =
			rho = calc_rhotp(x,rhol,rhog)
			mu = calc_mutp(x,mul,mug)
			phisq = calc_phi2(x,mul,mug,rhol,rhog)
		Poutk = calc_Pout(Pinj[j], mtot, tube, dii, x, xin, c, K, rhol,rhog,phisq,gamma,mul)) #Pacc will be zero for this
		#check syntax of do-while loop in python
		while True: 	
			Poutk1 = Poutk		
			TMLI = Tb + 0.5
			while True:
				calc_Qrad(Aio, Aoi, rio, roi, Twio, TMLI) + calc_Qspacer(spacer, dio, Twio, Twoi) 
				if (Qconv - Qrad - Qspacer < 10**(-12)): #checking equality
					break
			hl = 
			hg =
			Qtot = 
			mtot = ((1-x)*rhol + x*rhog)*c*calc_area(dii)
			hout = hin[j] + Qtot/mtot -0.5*(cout**2-c**2)- g*np.sin(gamma)*l
			xout = (hout - hl)/(hg - hl)
			Poutk =  calc_Pout(Pinj[j], mtot, tube, dii, xout, xin, c, K, rhol,rhog,phisq,gamma,mul))
		if (abs(Poutk - Poutk1) < 10**(-12)): #in MPa
			break
		j = j+1
		print(j)
		lz[j] = lz[j-1] + l
		hin[j] = hout
		Pinj[j] = Poutk
		TLHein = TLHeout[j-1]
		
		xin = x
	rhoout = rho	
	k = k+1
	coutn1 = mHetot/(rhoout*calc_area(dtubeout))) 
	ptv = 0.0975*np.exp(59.59*mRG) #mRG to be found
	n = n+1
	print(n)
	coutn = np.sqrt((psv + rhol*g*hLI + rhog*g*hLg - pout - rhoout*g*hsv - pfr - pfit)*2/rhoout)	#many unknowns in this formula to be found 	
	mHetot = coutn*rhoout * calc_area(dtubeout) 
	if (abs(coutn - coutn1)<0.01):
		break
#plotting values (p,h,rho,T)
plt.plot(lz,Pinj)
