from math import *
import math
from scipy.optimize import fsolve, minimize
import numpy as np
import pandas
import os
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Import file containing synthetic gases that match surface gas data
filename = 'Poas_model'
matchdatafilename = 'Poas_model_Match_favorite2gases.xlsx'
MatchData = pandas.read_excel(matchdatafilename, sheetname='All Data') #Read in gas data

##-------------SOME BOOKKEEPING-------------##
minP = 1.0
maxP = 10000.0
Pstep = 5.0
Prange = np.arange(minP, maxP, Pstep)

verbose = True #set to true to print out debugging statements in terminal
verbose_errors = True #set to true to print all calculated values thus far when an exception is raised

dQFMvalue = 0 #for model types 2 and 3, need to know fO2 relative to QFM buffer

MatchData["tempK"] = MatchData["temp"] + 273.15


##------------METHODS TO CLACULATE FUGACITIES--------------##
def Calc_fH2O_from_XH2O(thermo):
	fH2O = thermo["fH2Opure"] * thermo["XH2O"]

	return fH2O

def Calc_fH2(thermo):
	fH2 = thermo["fH2O"] / ( 10**Calc_KH2O(thermo["temp"]) * sqrt((10**(thermo["logfO2"]))) )
	return fH2

def Calc_fS2(thermo): #Equation 7 in Iacovino (2015) EPSL
	a2 = thermo["gammaS2"]**2
	b2 = thermo["KSO2"]**2
	c2 = thermo["fO2"]**2
	d2 = thermo["gammaSO2"]**2
	e2 = thermo["KH2S"]**2
	g2 = thermo["fH2"]**2
	h2 = thermo["gammaH2S"]**2

	v1 = 1.0 / (2 * d2 * h2)
	v2 = -(thermo["gammaS2"]**(1.5)) * (thermo["KSO2"] * thermo["fO2"] * thermo["gammaH2S"] + thermo["gammaSO2"] * thermo["KH2S"] * thermo["fH2"])

	s1 = thermo["gammaS2"] * b2 * c2 * h2
	s2 = 2 * thermo["gammaS2"] * thermo["KSO2"] * thermo["fO2"] * thermo["gammaSO2"] * thermo["KH2S"] * thermo["fH2"] * thermo["gammaH2S"]
	s3 = thermo["gammaS2"] * d2 * e2 * g2 + 4 * d2 * h2 * thermo["PStot"]
	sq = sqrt(s1 + s2 + s3)

	v31 = a2 * b2 * c2 * h2
	v32 = 2 * a2 * thermo["KSO2"] * thermo["fO2"] * thermo["gammaSO2"] * thermo["KH2S"] * thermo["fH2"] * thermo["gammaH2S"]
	v33 = a2 * d2 * e2 * g2
	v34 = 2 * thermo["gammaS2"] * d2 * h2 * thermo["PStot"]
	v3 	= v31 + v32 + v33 + v34

	fS2 = v1 * (v2 * sq + v3)

	return fS2

def Calc_fSO2(thermo):
	fSO2 = thermo["KSO2"] * sqrt(thermo["fS2"]) * thermo["fO2"]

	return fSO2

def Calc_fH2S(thermo):
	fH2S = thermo["KH2S"] * thermo["fH2"] * sqrt(thermo["fS2"])

	return fH2S

def Calc_fCO2(thermo):
	fCO2 = thermo["PCO2"] * thermo["gammaCO2"]

def Calc_fCO2_from_XCO2(thermo):
	fCO2 = thermo["fCO2pure"] * thermo["XCO2"]

	return fCO2

def Calc_fCO(thermo):
	fCO = thermo["fCO2"] / (thermo["KCO2"] * sqrt(thermo["fO2"]))

	return fCO

def Calc_All_Pure_fs(thermo):
	thermo["fCO2pure"] = thermo["press"] * thermo["gammaCO2"]
	thermo["fCOpure"] = thermo["press"] * thermo["gammaCO"]
	thermo["fSO2pure"] = thermo["press"] * thermo["gammaSO2"]
	thermo["fH2Spure"] = thermo["press"] * thermo["gammaH2S"]
	thermo["fH2Opure"] = thermo["press"] * thermo["gammaH2O"]

##-----------METHODS TO CALCULATE EQUILIBRIUM CONSTANTS-----------##
def Calc_KH2O(temp):
	# From Robie and Hemingway (1995). Text table
	# was fit with sixth-order polynomial to arrive
	# at the following coefficients.
	
	# H2 + 1/2O2 = H2O
	logK = 		(3.3426*10**(-17) 	* 	temp**6
			+	-2.40813*10**(-13)	*	temp**5
			+	7.10612*10**(-10)	*	temp**4
			+	-1.10671*10**(-06)	*	temp**3
			+	9.76038*10**(-04)	*	temp**2
			+	-4.84445*10**(-01)	*	temp
			+	1.21953*10**(2))
	return logK

def Calc_KH2S(temp):
	# From Robie and Hemingway (1995). Text table
	# was fit with sixth-order polynomial to arrive
	# at the following coefficients.
	
	# H2 + 1/2S2 = H2S
	logK = 		(1.882737*10**(-18) * temp**6
			+	-1.779266*10**(-14) * temp**5
			+	 6.319209*10**(-11) * temp**4
			+	-1.092048*10**(-07) * temp**3
			+	 9.824774*10**(-05) * temp**2
			+	-4.805344*10**(-02) * temp
			+	 1.389857*10**(01))
	return logK

def Calc_KSO2(temp):
	# From Robie and Hemingway (1995). Text table
	# was fit with sixth-order polynomial to arrive
	# at the following coefficients.
	
	# S2 + 1/2O2 = SO2
	logK =		 (4.01465*10**(-17) * temp**6
			+	-2.93845*10**(-13) 	* temp**5
			+	 8.78818*10**(-10) 	* temp**4
			+	-1.38066*10**(-06) 	* temp**3
			+	 1.21978*10**(-03) 	* temp**2
			+	-6.03476*10**(-01) 	* temp
			+	 1.54350*10**(02))
	return logK

def Calc_KCO2(temp):
	# From Wagman et al (1945). Text table
	# was fit with sixth-order polynomial to arrive
	# at the following coefficients.
	
	# CO + 1/2O2 = CO2
	logK =		 (9.11899*10**(-17)	* 	temp**6
			+	-5.65762*10**(-13) 	* 	temp**5
			+	 1.44706*10**(-9) 	*	temp**4
			+	-1.96925*10**(-06) 	*	temp**3
			+	 1.53277*10**(-03) 	*	temp**2
			+	-6.79138*10**(-01) 	*	temp
			+	 1.53288*10**(02))
	return logK

def Calc_all_logK(thermo):
	thermo["logK_H2O"] = Calc_KH2O(thermo["tempK"])
	thermo["logK_H2S"] = Calc_KH2S(thermo["tempK"])
	thermo["logK_SO2"] = Calc_KSO2(thermo["tempK"])
	thermo["logK_CO2"] = Calc_KCO2(thermo["tempK"])

def logK_to_K(thermo):
	thermo["KCO2"] = 10**thermo["logK_CO2"]
	thermo["KH2O"] = 10**thermo["logK_H2O"]
	thermo["KH2S"] = 10**thermo["logK_H2S"]
	thermo["KSO2"] = 10**thermo["logK_SO2"]

	if thermo["KCO2"] <= 0 or \
		thermo["KH2O"] <= 0 or \
		thermo["KH2S"] <= 0 or \
		thermo["KSO2"] <= 0:
		print "KCO2 = " + str(thermo["KCO2"])
		print "KH2O = " + str(thermo["KH2O"])
		print "KH2S = " + str(thermo["KH2S"])
		print "KSO2 = " + str(thermo["KSO2f"])
		if verbose_errors == True:
			print "\n"
			print "Verbose Errors:"
			print thermo
			print "\n"
		raise ValueError, "A calculated K value is negative or zero."
	else:
		pass


##-----------METHODS TO CALCULATE FUGACITY COEFFICIENTS-----------##
#Critical parameters cP, cT, o for relevant species
CPCO = {	"cT": 	133.15,
			"cP": 	34.9571,
			"o": 	0.049
		}

CPCO2 = {	"cT": 	304.15,
			"cP": 	73.8659,
			"o": 	0.225
		}

CPH2 = {	"cT": 	33.25,
			"cP": 	12.9696,
			"o": 	-0.218
		}

CPO2 = {	"cT": 	154.75,
			"cP": 	50.7638,
			"o": 	0.021
		}

CPS2 = {	"cT": 	208.15,
			"cP": 	72.954,
			"o": 	0.0 #need omega value for S2
		}

CPSO2 = {	"cT": 	430.95,
			"cP": 	778.7295,
			"o": 	0.0256
		}

CPH2S = {	"cT": 	373.55,
			"cP": 	90.0779,
			"o": 	0.081
		}

CPH2O = {	"cT": 	647.25,
			"cP": 	221.1925,
			"o": 	0.334
		}

##-----------METHODS TO CALCULATE PARTIAL PRESSURES-----------##

def Calc_PCO(thermo):
	PCO = thermo["fCO"] / thermo["gammaCO"]

	return PCO

def Calc_PCO2(thermo):
	PCO2 = thermo["fCO2"] / thermo["gammaCO2"]

	return PCO2

def Calc_PH2(thermo):
	PH2 = thermo["fH2"] / thermo["gammaH2"]

	return PH2

def Calc_PH2O(thermo):
	PH2O = thermo["fH2O"] / thermo["gammaH2O"]

	return PH2O

def Calc_PH2S(thermo):
	PH2S = thermo["fH2S"] / thermo["gammaH2S"]

	return PH2S

def Calc_PO2(thermo):
	PO2 = thermo["fO2"] / thermo["gammaO2"]

	return PO2

def Calc_PS2(thermo):
	PS2 = thermo["fS2"] / thermo["gammaS2"]

	return PS2

def Calc_PSO2(thermo):
	PSO2 = thermo["fSO2"] / thermo["gammaSO2"]

	return PSO2

def Calc_PStot(thermo):
	PStot = thermo["press"] - thermo["PCO"] - thermo["PCO2"] - thermo["PH2"] - thermo["PH2O"] - thermo["PO2"]

	return PStot

def Calc_all_PartialPressures(thermo):
	thermo["PH2"] = Calc_PH2(thermo)
	thermo["PH2O"] = Calc_PH2O(thermo)
	thermo["PH2S"] = Calc_PH2S(thermo)
	thermo["PO2"] = Calc_PO2(thermo)
	thermo["PS2"] = Calc_PS2(thermo)
	thermo["PSO2"] = Calc_PSO2(thermo)

	thermo["PCtot"] = 	(thermo["press"] 
							- thermo["PH2"] 
							- thermo["PH2O"] 
							- thermo["PH2S"]
							- thermo["PO2"]
							- thermo["PS2"]
							- thermo["PSO2"]
						)
											
	thermo["PCO"] = thermo["gammaCO2"] * thermo["PCtot"] / (thermo["gammaCO2"] + thermo["gammaCO"] * thermo["KCO2"] * sqrt(thermo["fO2"]))

	thermo["PCO2"] = (thermo["press"]
						- thermo["PCO"]
						- thermo["PH2"] 
						- thermo["PH2O"] 
						- thermo["PH2S"]
						- thermo["PO2"]
						- thermo["PS2"]
						- thermo["PSO2"]
					)

##-----------METHODS TO CALCULATE FLUID MOL FRACTIONS-----------##

def Calc_XCO(thermo):
	XCO = thermo["fCO"] / (thermo["gammaCO"] * thermo["press"])

	return XCO

def Calc_XCO2(thermo): #ONLY USE if XCO2 is not given as a user input
	XCO2 = thermo["fCO2"] / (thermo["gammaCO2"] * thermo["press"])

	return XCO2

def Calc_XH2O(thermo):
	XH2O = thermo["fH2O"] / (thermo["gammaH2O"] * thermo["press"])

	return XH2O

def Calc_XH2(thermo):
	XH2 = thermo["fH2"] / (thermo["gammaH2"] * thermo["press"])

	return XH2

def Calc_XO2(thermo):
	XO2 = thermo["fO2"] / (thermo["gammaO2"] * thermo["press"])

	return XO2

def Calc_XS2(thermo):
	XS2 = thermo["fS2"] / (thermo["gammaS2"] * thermo["press"])

	return XS2

def Check_XCO2_value(thermo): #Checks that user input XCO2 value is not too large
	max_val = 1 - thermo["XCO"] - thermo["XH2O"] - thermo["XH2"] - thermo["XO2"]

	while thermo["XCO2"] >= max_val:
		thermo["XCO2"] = max_val
		thermo["fCO2"] = Calc_fCO2_from_XCO2(thermo) #recalc fCO2 using new XCO2
		thermo["fCO"] = Calc_fCO(thermo) #recalc fCO using new fCO2
		thermo["XCO"] = Calc_XCO(thermo) #recalc XCO using new fCO
		max_val = 1 - thermo["XCO"] - thermo["XH2O"] - thermo["XH2"]
		print "Warning: User input XCO2 value is too large. Using newly calulated value of " + str(max_val)



##-----------METHODS TO CALCULATE FUGACITY COEFFICIENTS-----------##
#Redlich Kwong Equation of State
#Code derived from http://people.ds.cam.ac.uk/pjb10/thermo/pure.html - Patrick J. Barrie 30 October 2003
def RedlichKwong(temp, press, CP): #takes temperature in C and converts to K; takes pressure in bar and uses bar; CP = critical parameters in a dict as cP, cT, omega
	tempK = temp + 273.15
	R = 8.3145

	#Calculate a and b parameters (depend only on critical parameters)...
	a = 0.42748 * R**2 * CP["cT"]**(2.5) / (CP["cP"] * 10**5)
	b = 0.08664 * R * CP["cT"] / (CP["cP"] * 10**5)
	kappa = 0.0

	#Calculate coefficients in the cubic equation of state...
	#coeffs: (C0, C1, C2, A, B)
	A = a * press * 10**5 / (sqrt(tempK) * (R * tempK)**2)
	B = b * press * 10**5 / (R * tempK)
	C2 = -1.0
	C1 = A - B - B * B
	C0 = -A * B

	#Solve the cubic equation for Z0 - Z2, D...
	Q1 = C2 * C1 / 6 - C0 / 2 - C2**3 / 27
	P1 = C2**2 / 9 - C1 / 3
	D = Q1**2 - P1**3

	if D >= 0:
		kOneThird = 1.0 / 3.0

		absQ1PSqrtD = fabs(Q1 + sqrt(D))
		temp1 = absQ1PSqrtD**kOneThird
		temp1 *= (Q1 + sqrt(D)) / absQ1PSqrtD

		absQ1MSqrtD = fabs(Q1 - sqrt(D))
		temp2 = absQ1MSqrtD**kOneThird
		temp2 *= (Q1 - sqrt(D)) / absQ1MSqrtD

		Z0 = temp1 + temp2 - C2 / 3
	else:
		temp1 = Q1**2 / (P1**3)
		temp2 = sqrt(1.0 - temp1) / sqrt(temp1)
		temp2 *= Q1 / fabs(Q1)

		gamma = atan(temp2)

		if gamma < 0:
			gamma = gamma + pi 

		Z0 = 2 * sqrt(P1) * cos(gamma/3) - C2 / 3
		Z1 = 2 * sqrt(P1) * cos((gamma + 2 * pi) / 3) - C2/3
		Z2 = 2 * sqrt(P1) * cos((gamma + 4 * pi) / 3) - C2/3

		if Z0 < Z1:
			temp0 = Z0
			Z0 = Z1
			Z1 = temp0

		if Z1 < Z2:
			temp0 = Z1
			Z1 = Z2
			Z2 = temp0

		if Z0 < Z1:
			temp0 = Z0
			Z0 = Z1
			Z1 = temp0

	#Determine the fugacity coefficient of first root and departure functions...
	#calcdepfns(coeffs[3], 	coeffs[4], 	paramsab[0], 	Z[0])
	#calcdepfns(A, 			B, 			kappa, 			Z)

	#Calculate Departure Functions
	gamma = exp(Z0 - 1 - log(Z0-B) - A * log(1+B/Z0)/B)

	Hdep = R * tempK * (Z0 - 1 - 1.5*A*log(1+B/Z0)/B)

	Sdep = R * (log(Z0-B) - 0.5*A*log(1+B/Z0)/B)

	return gamma

def CalcAllGammas(thermo):
	for i in range(len(fugacity_coefficients_list)):
		thermo[fugacity_coefficients_list[i]] = RedlichKwong(thermo["temp"], thermo["press"], CP_list[i])
		if thermo[fugacity_coefficients_list[i]] <= 0:
			print thermo[fugacity_coefficients_list[i]]
			if verbose_errors == True:
				print "\n"
				print "Verbose Errors:"
				print thermo
			else:
				pass
			raise ValueError, "The fugacity coefficient (gamma value) for " + str(thermo[fugacity_coefficients_list[i]]) + "is negative or zero."
		else:
			pass


##-----------CALC FO2 FROM DELTA QFM VALUE------------##

def fO2_from_dQFM(thermo):
	fO2_at_QFM = -25096.3/thermo["tempK"] + 8.735 + 0.11*(thermo["press"] - 1)/thermo["tempK"]
	new_log_fO2 = fO2_at_QFM + thermo["dQFMvalue"]
	new_fO2 = 10**(new_log_fO2)

	return new_fO2, new_log_fO2

##------------RESPECIATE FLUIDS AT LOW P--------------##

def respeciate(thermo, Pressure):
	XHtot = thermo["XH2Otot"] * 2 / 3
	XStot = thermo["XStot"]
	XCtot = thermo["XCO2tot"] * 1 / 3
	XOtot = thermo["XH2Otot"] * 1 / 3 + thermo["XCO2tot"] * 2 / 3

	B = thermo["gammaH2"]
	P = Pressure
	C = thermo["KH2O"]
	D = thermo["fO2low"]
	sD = sqrt(D)
	E = thermo["gammaH2O"]
	F = thermo["KH2S"]
	G = thermo["gammaH2S"]
	J = thermo["gammaS2"]
	K = thermo["KSO2"]
	L = thermo["gammaSO2"]
	M = thermo["KCO2"]
	N = thermo["gammaCO2"]
	Q = thermo["gammaCO"]

	#FIRST calculate fH2 and fS2 using fsolve, two eqns; two unknowns (eqn 9 in Iacovino, 2015)
	def equations(p):
		fH2, fS2 = p
		return 	(
					( (fH2/(B*P)) 	+ ((2 * C * fH2 * sD)/(3 * E * P))		+ ((2 * F * fH2 * sqrt(abs(fS2)))/(3 * G * P)) 	- XHtot), 
					( (fS2/(J * P)) 	+ ((F * fH2 * sqrt(abs(fS2)))/(3 * G * P))	+ ((K * D * sqrt(abs(fS2)))/(3 * L * P))		- XStot)
				)

	fH2, fS2 = fsolve(equations, (0.00001, 0.00001))

	if fS2 <= 0 or fH2 <=0:
		fS2calcd = fS2
		fH2calcd = fH2
		fS2 = 0.001
		fH2 = 0.001
		if verbose_errors == True:
			print "Warning: Calculated negative fS2 or fH2. Not saving this result. Moving on to next calcualtion..."
			print "fS2: " + str(fS2calcd)
			print "fH2: " + str(fH2calcd)
		else:
			pass
		thermo["dont_save_this_iteration"] = True
	else:
		pass
	#SECOND calculate fCO (eqn 10 in Iacovino, 2015)
	def fCO_func(fCO):
		return (((M * fCO * sD)/(3 * N * Pressure)) + ((fCO)/(2 * Q * Pressure))	- XCtot)

	[fCO] = fsolve(fCO_func, 0.001)

	#THIRD calculate fCO2 using calc'd fCO and known fO2 value
	fCO2 = M * fCO * sD

	#FOURTH calcualte fSO2 using calc'd fS2 and known fO2 value
	fSO2 = K * sqrt(fS2) * D

	#FIFTH calculate fH2S using calc'd fH2 and fS2 values
	fH2S = F * fH2 * sqrt(fS2)

	#SIXTH calculate fH2O using calc'd fH2 and knwn fO2 value
	fH2O = C * sD * fH2

	thermo["fH2"] 	= fH2
	thermo["fS2"] 	= fS2
	thermo["fCO"] 	= fCO
	thermo["fCO2"] 	= fCO2
	thermo["fSO2"]	= fSO2
	thermo["fH2S"] 	= fH2S
	thermo["fH2O"]	= fH2O

	#sync all fO2 values to only use input Pressure value moving forward
	thermo["fO2"] = D

	if 	thermo["fH2"] <= 0 or \
		thermo["fS2"] <= 0 or \
		thermo["fCO"] <= 0 or \
		thermo["fCO2"] <= 0 or \
		thermo["fSO2"] <= 0 or \
		thermo["fH2S"] <= 0 or \
		thermo["fH2O"] <= 0:

		print "XCtot = " + str(XCtot)
		print "fH2 = " + str(thermo["fH2"])
		print "fS2 = " + str(thermo["fS2"])
		print "fCO = " + str(thermo["fCO"])
		print "fCO2 = " + str(thermo["fCO2"])
		print "fSO2 = " + str(thermo["fSO2"])
		print "fH2S = " + str(thermo["fH2S"])
		print "fH2O = " + str(thermo["fH2O"])
		if verbose_errors == True:
			print "\n"
			print "Verbose Errors:"
			print thermo
		else:
			pass
		raise ValueError, "A fugacity is negative or zero."
	else:
		Calc_all_Xs_from_fS(thermo)

##------------CALCULATE MOLE FRACS FROM FUGACITIES AT NEW PRESSURE--------------##

def Calc_all_Xs_from_fS(thermo):
	thermo["XCO"] 	= Calc_XCO(thermo)
	thermo["XCO2"] 	= Calc_XCO2(thermo)
	thermo["XH2O"] 	= Calc_XH2O(thermo)
	thermo["XH2"]	= Calc_XH2(thermo)
	thermo["XH2S"]	= thermo["fH2S"] / (thermo["gammaH2S"] * thermo["press"])
	thermo["XO2"]	= Calc_XO2(thermo)
	thermo["XS2"]	= Calc_XS2(thermo)
	thermo["XSO2"]	= thermo["fSO2"] / (thermo["gammaSO2"] * thermo["press"])

	if 	thermo["XH2"] <= 0 or \
		thermo["XS2"] <= 0 or \
		thermo["XCO"] <= 0 or \
		thermo["XCO2"] <= 0 or \
		thermo["XSO2"] <= 0 or \
		thermo["XH2S"] <= 0 or \
		thermo["XH2O"] <= 0:
		print "XH2 = " + str(thermo["XH2"])
		print "XS2 = " + str(thermo["XS2"])
		print "XCO = " + str(thermo["XCO"])
		print "XCO2 = " + str(thermo["XCO2"])
		print "XSO2 = " + str(thermo["XSO2"])
		print "XH2S = " + str(thermo["XH2S"])
		print "XH2O = " + str(thermo["XH2O"])
		if verbose_errors == True:
			print "\n"
			print "Verbose Errors:"
			print thermo
		else:
			pass
		raise ValueError, "A respeciated fluid mole fraction is negative or zero."
	else:
		pass

##------------NORMALIZE FINAL X VALUES--------------##
def normalize_final_Xs(thermo):
	Normed_Xs = {}
	orig_X_sum = (	thermo["XCO"] + 
					thermo["XCO2"] +
					thermo["XH2"] + 
					thermo["XH2O"] + 
					thermo["XH2S"] +
					thermo["XO2"] +
					thermo["XS2"] +
					thermo["XSO2"])

	Normed_Xs["XCO"] 	= 	thermo["XCO"] / orig_X_sum
	Normed_Xs["XCO2"]	=	thermo["XCO2"] / orig_X_sum
	Normed_Xs["XH2"]	=	thermo["XH2"] / orig_X_sum
	Normed_Xs["XH2O"]	=	thermo["XH2O"] / orig_X_sum
	Normed_Xs["XH2S"]	=	thermo["XH2S"] / orig_X_sum
	Normed_Xs["XO2"]	=	thermo["XO2"] / orig_X_sum
	Normed_Xs["XS2"]	=	thermo["XS2"] / orig_X_sum
	Normed_Xs["XSO2"]	=	thermo["XSO2"] / orig_X_sum

	Normed_Xs["X_Sum"] = sum(Normed_Xs.values())

	for key in Normed_Xs:
		thermo[key] = Normed_Xs[key]

	return Normed_Xs

##-------------DEBUGGING THE HARD WAY-------------##
def print_all_Xs(thermo):
	print "XCO = " + str(thermo["XCO"])
	print "XCO2 = " + str(thermo["XCO2"])
	print "XH2 = " + str(thermo["XH2"])
	print "XH2O = " + str(thermo["XH2O"])
	print "XH2S = " + str(thermo["XH2S"])
	print "XO2 = " + str(thermo["XO2"])
	print "XS2 = " + str(thermo["XS2"])
	print "XSO2 = " + str(thermo["XSO2"])
	X_Sum = (thermo["XCO"]
				+ thermo["XCO2"]
				+ thermo["XH2"]
				+ thermo["XH2O"]
				+ thermo["XH2S"]
				+ thermo["XO2"]
				+ thermo["XS2"]
				+ thermo["XSO2"])
	thermo["X_Sum"] = X_Sum
	print "X sum = " + str(X_Sum)

##-------------Some Global Variables----------------##
#Calculate fugacity coefficients and put them in the thermo dict
fugacity_coefficients_list = ['gammaCO', 'gammaCO2', 'gammaH2', 'gammaH2O', 'gammaH2S', 'gammaO2', 'gammaS2', 'gammaSO2'] #List of all fugacity coefficients
CP_list = [CPCO, CPCO2, CPH2, CPH2O, CPH2S, CPO2, CPS2, CPSO2] #List of all critical parameter dicts

##--------------------------DEPRESSURIZE-----------------------##
#Degass at 1 bar steps
depressurizelist = []
depressurize_new_fluids_list = []
gasnumber = 1
for l in Prange:
	for index, row in MatchData.iterrows():
		#SET PRESSURE FOR THIS LOOP
		row["press"] = l

		row["Gas_ID"] = gasnumber
		row["dQFMvalue"] = dQFMvalue
		row["fO2low"], row["logfO2"] = fO2_from_dQFM(row)
		print gasnumber

		#run the model.....
		row["dont_save_this_iteration"] = False
		Calc_all_logK(row)
		CalcAllGammas(row)
		logK_to_K(row)
		respeciate(row, l)
		

		Calc_all_Xs_from_fS(row)
		Normd_fluid = normalize_final_Xs(row)

		if row["dont_save_this_iteration"] == False:
			Normd_fluid["ratio_SO2_CO2"] = Normd_fluid["XSO2"] / Normd_fluid["XCO2"]
			Normd_fluid["ratio_H2S_SO2"] = Normd_fluid["XH2S"] / Normd_fluid["XSO2"]

			#save data externally to this loop
			depressurize_new_fluids_list.append(Normd_fluid)
			row["ID_number"] = gasnumber
			depressurizelist.append(row) 

	gasnumber += 1

DepressurizeFluidFrame = pandas.DataFrame(depressurizelist)

writer = pandas.ExcelWriter(filename + '_depressurize-steps.xlsx', engine='xlsxwriter')
DepressurizeFluidFrame.to_excel(writer)
writer.save() #Close the Pandas Excel writer and output the Excel file

print "Success! Saved file " + str(filename) + "_depressurize-steps.xlsx"

##==============================================================##