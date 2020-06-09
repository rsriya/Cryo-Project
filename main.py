#To find
#1. Subroutine to calculate Twio
#2. How to calculate Pin
#3. How to find rho and mu


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
do:
#check syntax of do-while loop in python
	j = 0
	Pinj[0] = Pin
	while (j<integer(ltot/l)):
		K = 0 #no consideration of fitting loss
		#calculate K by the following if conditions if want to consider		
		if lz< :
			tube = s
			dii =
		elif lz<:
			tube = r
			dii =
		else:
			tube = s
			dii =
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
		Re = rho * cin * dii/mu #check formula
		Poutk = calc_Pout(Pinj[j], mtot, tube, dii, Re, x, xin, c, K, rhol,rhog,phisq)) #Pacc will be zero for this
		#check syntax of do-while loop in python
		do: #in MPa	
			Poutk1 = Poutk		
			TMLI = Tb + 0.5
			hl = 
			hg =
			
			x = (hout - hl)/(hg - hl)
			Poutk =  calc_Pout(Pinj[j], mtot, tube, dii, Re, x, xin, c, K, rhol,rhog,phisq))
		}while (abs(Poutk - Poutk1) > 10**(-12))
		j = j+1
		print(j)
		lz[j] = lz[j-1] + l
		hin = hout[j-1]
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
}while (abs(coutn - coutn1)>0.01)
#plotting values (p,h,rho,T)
plt.plot(lz,Pinj)
