from constants import * 
import math
import numpy as np

def calc_QconvA(NuLO, Twii, dii):
	QconvA = NuLO * lambdal * (Twii-TB)/dii
	return QconvA

def calc_QBA(cpl,Twii,Tb,Prl,r):
	QBA = ( cpl*(Twii-Tb)/(Prl**(1.7)*r*0.169) )**3 / np.sqrt(sigmal/(g*(rhol-rhog)) ) * mul * r
	return QBA

def calc_QSNBA(NuLO, Twii, dii, cpl,Twii,Tb,Prl,r):
	QSNBA = calc_QconvA(NuLO, Twii, dii) + calc_QBA(cpl,Twii,Tb,Prl,r)
	return QSNBA

def calc_phi(rio, roi):
	y1 = rio/roi
	y2 = l/roi 
	zeta = y2**2 - y1**2 + 1
	sai = y2**2 + y1**2 - 1
	phi = 1/y1 - 1/(np.pi * y1) * (math.acos(zeta/sai) - 1/(2*y2) * (np.sqrt((sai + 2)**2 - 4*y1**2) * math.acos(zeta/(sai*y1)) + zeta*math.asin(1/y1) - np.pi/2*sai ))
	return phi

def calc_Rtube(rii, roi, h):
	Rtube = 1/(2*np.pi*l) * (1/(h*rii) + 1/lambdatube * np.log(roi/rii))

def calc_Xtt(x):
	Xtt = (rhog/rhol)**0.5 * (mul/mug)**0.1 * ((1-x)/x)**0.9
	return Xtt

def calc_hsp(Re,Pr,lamda,dii):
	hsp = 0.023 * Re**0.8 ** Pr**(1/3) * lamda/dii
	return hsp

def calc_hFCE(x, NuLO, dii):
	hFCE = 3.5 * (calc_Xtt(x))**(-0.5) * NuLO * lambdal/dii
	return hFCE

def calc_hSNB(Twii):
	hSNB = calc_QSNBA(NuLO, Twii, dii, cpl,Twii,Tb,Prl,r) * (Twii - TB)**(-1) 
	return hSNB

def calc_Qrad(Aio, Aoi, rio, roi, Twio, TMLI):
	C = 5.67/(1/epsio - 1 + 1/calc_phi(rio, roi) + (1/epsoi-1)*Aio/Aoi)
	Qrad = C*Aio*((Twio/100)**4 - (TMLI/100)**4)
	return Qrad

def calc_QMLIrad(AMLIo, Tm, TMLI, Twoi):
	QMLI = AMLIo * (5.754 * 10**(-7) * Nbar**2.02 * Tm * (TMLI-Twoi)/(Ns+1) + 1.089*10**(-9) * epsTRAlu*epsTRMylar * (TMLI**4.67-Twoi**4.67)/((epsTRAlu+epsTRMylar)*Ns) )
	return QMLI 

def calc_QMLIcond(doi, TMLI, Twoi):
	QMLI = lambdaAlu*deltaAlu*l/(Ns*np.pi*(doi+deltaMLI)) * (TMLI - Twoi)
	return QMLI

def calc_Qtube(Twoi, rii, roi, h):
	Qtube = (Twoi - TB)/calc_Rtube(rii, roi, h)	
	return Qtube

def calc_QONB(A,deltaT,NuLO,r,dii,x):
	QONB = A*98*(rhol-rhog)*sigmal*lambdal*deltaT*NuLO**2/(rhol*rhog*r*calc_Xtt(x)*dii**2)
	return QONB

def calc_Qspacer(dio, Twio, Twoi):
	Qspacer = (dio/(2*lambdaPTFE)*np.log(dio/(doi+2*deltaMLI)) )**(-1) * np.pi*dio*deltaPTFE*0.1*(Twio - Twoi)
	return Qspacer

