from constants import * 
import numpy as np
from scipy import optimize

def calc_rhotp(x,rhol,rhog):
	return 1/(x/rhog + (1-x)/rhol)

def calc_mutp(x,mul,mug):
	return 1/(x/mug + (1-x)/mul)

def calc_alpha(x,rhol,rhog):
	return 1/(1 + (1-x)/x * rhog/rhol)

def calc_f_rigid(dii, Re):
	#a function in terms of f whose root has to be found
	return 1 + 2 * (f ** 0.5) * np.log10( ktube/(3.7*dii) + 2.51/(Re*f**0.5) )

def calc_f_flex(dii, Re):
	a = 0.01588 * (s-e)/t - 0.00215
	b = 0.2987 * t*e/s**2 - 0.0313
	f = (1 + 7.898 * dii/rbend) * 4 * a * Re ** b #(7) 


def calc_Pfr(tube, dii, x, c, rhol, rhog, mul, phisq):
	G = ((1-x)*rhol + x*rhog)*c
	RLO = G*dii/mul
	if tube == 'r':
		f = optimize.root_scalar(calc_f_rigid, args = (dii, ReLO))
	else:
		f = calc_f_flex(dii, ReLO)
	Pfr = (phisq * f * G ** 2 * l) / (2 * dii * rhol)
	return Pfr

def calc_Pg(x,rhol,rhog, gamma): #Pressure loss
	return calc_rhotp(x,rhol,rhog)*g*np.sin(gamma)*l 

def calc_Pacc(G, xout, xin, rhol,rhog):  
	alphain = calc_alpha(xin,rhol,rhog)
	alphaout = calc_alpha(xout,rhol,rhog)
	win = ((1-xin)**2)/((1-alphain)*rhol) + (xin**2)/(alphain*rhog)
	wout = ((1-xout)**2)/((1-alphaout)*rhol) + (xout**2)/(alphaout*rhog)
	return G**2 * (wout - win)

def calc_Pfitt(c,K,rhol,phisq):
	return phisq**0.5 * K * 0.5 * rhol * c**2

def calc_Pout(Pin, tube, dii, xout, xin, c, K, rhol,rhog,phisq, gamma,mul):
	return Pin - (calc_Pfr(tube, dii, xin, c, rhol, rhog, mul, phisq) + calc_Pg(xin,rhol,rhog, gamma) + calc_Pacc(G, xout, xin, rhol, rhog) + calc_Pfitt(c,K,rhol,phisq))


 

