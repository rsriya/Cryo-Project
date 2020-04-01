from constants import * 
import numpy as np
from scipy import optimize

def calc_rhotp(x):
	return 1/(x/rhog + (1-x)/rhol)

def calc_mutp(x):
	return 1/(x/mug + (1-x)/mul)

def calc_alpha(x):
	return 1/(1 + (1-x)/x * rhog/rhol)

def calc_f_rigid(dii, Re):
	#a function in terms of f whose root has to be found
	return 1 + 2 * (f ** 0.5) * np.log10( ktube/(3.7*dii) + 2.51/(Re*f**0.5) )

def calc_f_flex(dii, Re):
	a = 0.01588 * (s-e)/t - 0.00215
	b = 0.2987 * t*e/s**2 - 0.0313
	f = (1 + 7.898 * dii/rbend) * 4 * a * Re ** b #(7) 

def calc_phi2(x):
	return (1 + x*(mul - mug)/mug)**(-0.25) * (1 + x*(rhol/rhog - 1))

def calc_Pfr(mtot, tube, dii, Re, x, G):
	if tube == 'r':
		f = optimize.root_scalar(calc_f_rigid, args = (dii, Re))
	else:
		f = calc_f_flex(dii, Re)
	Pfr = (calc_phi2(x) * f * G ** 2 * l) / (2 * dii * rhol)
	return Pfr

def calc_Pg(x):
	return calc_rhotp(x)*g*np.sin(gamma)*l 

def calc_Pacc(G, x, xout): #assuming xin = x
	alphain = calc_alpha(x)
	alphaout = calc_alpha(xout)
	win = ((1-x)**2)/((1-alphain)*rhol) + (x**2)/(alphain*rhog)
	wout = ((1-xout)**2)/((1-alphaout)*rhol) + (xout**2)/(alphaout*rhog)
	return G**2 * (wout - win)

def calc_Pfitt(c,K):
	return calc_phi2(x)**0.5 * K * 0.5 * rhol * c**2

def calc_Pout(Pin, mtot, tube, dii, Re, x, xout, fit, c, K):
	Aii = np.pi * dii**2 / 4
	G = mtot/Aii
	return Pin - (calc_Pfr(mtot, tube, dii, Re, x, G) + calc_Pg(x) + calc_Pacc(G, x, xout) + calc_Pfitt(c,K))


 

