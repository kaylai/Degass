from math import *
import math
from scipy.optimize import fsolve, minimize
import numpy as np
import pandas
import os
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##-----------------USER INPUTS---------------##
#MELT INCLUSION WT% AS (H2O, CO2, S)
#Temp in degrees celsius
#Press in bars

filename = 'Poas_model'
save_plots = False
verbose = False #set to true to print out debugging statements in terminal
verbose_errors = True #set to true to print all calculated values thus far when an exception is raised
MatchBoth = True #set to true to only find matches that satisfy BOTH SO2/CO2 and H2S/SO2 ratios. False results in matches calculated for each ratio separately.

###For algorithm to match with surface gas data
gasdatafilename = "PreppedPoasData.xlsx" #filename of excel file with surface gas data
tolerance = 0.20 #percent allowable tolerance around calcaulted values to find a match with observed data. 0.05 = 5%

dQFMvalue = 0 #for model types 2 and 3, need to know fO2 relative to QFM buffer

##############user inputs for matching model
#####values as dissolved wt% concentrations in deep melt inclusions
minH2O = 1.0
maxH2O = 6.0
H2Ostep = 0.5
minCO2 = 0.03
maxCO2 = 0.30
CO2step = 0.01
minS = 0.03
maxS = 0.5
Sstep = 0.05
minP = 1.0
maxP = 2000.0
Pstep = 1.0

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

##-------------------SOME BOOK KEEPING---------------##

###Some common assertions
assert minH2O > 0, "Minimum dissolved H2O in deep melt inclusion cannot be zero. Try 0.0001 if you want a very low value."
assert minCO2 > 0, "Minimum dissolved CO2 in deep melt inclusion cannot be zero. Try 0.0001 if you want a very low value."
assert minS > 0, "Minimum dissolved S in deep melt inclusion cannot be zero. Try 0.0001 if you want a very low value."

###For algorithm to match with surface gas data
GasData = pandas.read_excel(gasdatafilename) #Read in gas data
minvalgasdata1 = GasData["SO2/CO2"].min()
maxvalgasdata1 = GasData["SO2/CO2"].max()
minvalgasdata2 = GasData["H2S/SO2"].min()
maxvalgasdata2 = GasData["H2S/SO2"].max()

##------------------CONSTANTS---------------##

MW = { "SiO2": 60.08,
		"TiO2": 79.866,
		"Al2O3": 101.96,
		"FeOstar": 71.844,
		"FeO": 71.844,
		"Fe2O3": 159.69,
		"MnO": 70.9374,
		"MgO": 40.3044,
		"CaO": 56.0774,
		"Na2O": 61.9789,
		"K2O": 94.2,
		"P2O5": 283.886,
		"H2O": 18.02,
		"CO2": 44.01,
		"S": 32.065,
		"O": 15.9994,
		"H": 1.00794,
		"C": 12.0107
	}


##----------------DEFINE METHODS----------------##
def Convert_wt_to_molfrac(meltcomp):
	MolPropOx = {}
	MolFrac = {}
	for key in meltcomp:
		MolPropOx[key] = meltcomp[key]/MW[key]

	MPO_sum = sum(MolPropOx.values())

	for key in MolPropOx:
		MolFrac[key] = (MolPropOx[key]/MPO_sum)

	#calculate volatiles separately for fluid phase composition
	H2O = MolFrac["H2O"]
	CO2 = MolFrac["CO2"]
	S = MolFrac["S"]
	sumfl = H2O + CO2 + S

	if sumfl <= 0:
		print "sumfl = " + str(sumfl)
		if verbose_errors == True:
			print "\n"
			print "Verbose Errors:"
			print MeltInclusion_deep_thermo
		else:
			pass
		raise ValueError, "Sume of fluid mole fractions is negative or zero."

	else:
		H2Ofl = H2O/sumfl
		CO2fl = CO2/sumfl
		Sfl = S/sumfl

		MolFrac["XH2Ofl"] = H2Ofl
		MolFrac["XCO2fl"] = CO2fl
		MolFrac["XSfl"] = Sfl

		return MolFrac

def Difference_MI_volatiles(deep_MI, shallow_MI):
	difference = {}
	for key in deep_MI:
		difference[key] = deep_MI[key]
		if key == "H2O" or key == "CO2" or key == "S":
			difference[key] = (deep_MI[key] - shallow_MI[key])
			
			#Fluid volatile values cannot be <= 0 or resulting fugacities will be <=0. This is partly protected with the assert that deep dissolved \
			#volatile values cannot be <=0, but if shallow dissolved > deep dissolved, the difference can be <=0. This protects against that by simply \
			#assigning it a new value of 0.0001 to simulate a very small value. Probably good to find a good way to notify the user this is happening.
			if difference[key] <= 0:
				difference[key] = 0.001
				MeltInclusion_deep_thermo["dont_save_this_iteration"] = True
			else:
				pass
		else:
			pass
	return difference


##------------METHODS TO CLACULATE FUGACITIES--------------##
def Calc_fH2O_Moore1998(meltcompx, thermo): #must pass melt comp in hydrous mole fraction, totalling 1 #thermo is the thermo data (P, T, fO2) for your sample
	a = 2656
	bAl2O3 = -1.997
	bFeOt = -0.9275
	bNa2O = 2.736
	c = 1.171
	d = -14.21

	press = thermo["press"]
	temp = thermo["temp"]

	XAl2O3 = meltcompx["Al2O3"]
	XFeOt = meltcompx["FeO"] + 0.8998 * meltcompx["Fe2O3"] + meltcompx["FeOstar"]
	XNa2O = meltcompx["Na2O"]
	XH2Om = meltcompx["H2O"]
	LNX = 2 * log(XH2Om)

	fugacity_H2O = exp( (LNX - (a/temp) - ((bAl2O3 * XAl2O3 * (press/temp)) + (bFeOt * XFeOt * (press/temp)) + (bNa2O * XNa2O * (press/temp))) - d)/c) #Moore1998 model

	return fugacity_H2O

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



##-----------METHODS TO CALCULATE FLUID RATIOS-----------##
def Calc_all_fluid_ratios(thermo, comp):
	thermo["fC_ratio"] = thermo["fCO"] / thermo["fCO2"]
	thermo["fCpure_ratio"] = thermo["fCOpure"] / thermo["fCO2"]
	thermo["XC_ratio"] = thermo["fC_ratio"] / thermo["fCpure_ratio"]

	thermo["fS_ratio"] = thermo["fH2S"] / thermo["fSO2"]
	thermo["fSpure_ratio"] = thermo["fH2Spure"] / thermo["fSO2pure"]
	thermo["XS_ratio"] = thermo["fS_ratio"] / thermo["fSpure_ratio"]

	XSO2_XStot = 1/(1 + thermo["XS_ratio"])
	XCO2_Ctot = 1/(1 + thermo["XC_ratio"])

	S_in_SO2 = XSO2_XStot * comp["S"]
	S_in_H2S = comp["S"] - S_in_SO2

	CO2_in_CO2tot = XCO2_Ctot * thermo["XCO2"]
	CO_in_CO2tot = thermo["XCO2"] - CO2_in_CO2tot

	thermo["XSO2"] = S_in_SO2 * ((MW["S"] + MW["O"] + MW["O"])/MW["S"])
	thermo["XH2S"] = S_in_H2S * ((MW["H"] + MW["H"] + MW["S"])/MW["S"])
	thermo["XCO"] = CO_in_CO2tot * ((MW["C"] + MW["O"])/(MW["C"] + MW["O"] + MW["O"]))
	thermo["XCO2"] = CO2_in_CO2tot



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

##-----------CALC FO2 FROM DELTA QFM VALUE------------##

def fO2_from_dQFM(thermo):
	fO2_at_QFM = -25096.3/thermo["tempK"] + 8.735 + 0.11*(thermo["press"] - 1)/thermo["tempK"]
	new_log_fO2 = fO2_at_QFM + thermo["dQFMvalue"]
	new_fO2 = 10**(new_log_fO2)

	return new_fO2, new_log_fO2

##------------RESPECIATE FLUIDS AT LOW P--------------##

def respeciate(thermo, lowP):
	XHtot = thermo["XH2Otot"] * 2 / 3
	XStot = thermo["XStot"]
	XCtot = thermo["XCO2tot"] * 1 / 3
	XOtot = thermo["XH2Otot"] * 1 / 3 + thermo["XCO2tot"] * 2 / 3

	B = thermo["gammaH2"]
	P = lowP
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
		fS2 = 0.001
		fH2 = 0.001
		if verbose_errors == True:
			print "Warning: Calculated negative fS2 or fH2. Not saving this result. Moving on to next calcualtion..."
		else:
			pass
		thermo["dont_save_this_iteration"] = True
	else:
		pass
	#SECOND calculate fCO (eqn 10 in Iacovino, 2015)
	def fCO_func(fCO):
		return (((M * fCO * sD)/(3 * N * lowP)) + ((fCO)/(2 * Q * lowP))	- XCtot)

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

	#sync all fO2 values to only use lowP value moving forward
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

##-------------TheModel----------------##
#For differencing melt inclusions and printing normalized mole fraction values

def TheModel(thermo):
	#reset solvefS2 debugger:
	thermo["solvefS2"] = 1
	#Cacltulate fO2 at low P from buffer
	#first, store deltaQFM value in thermodynamic dict
	thermo["dQFMvalue"] = dQFMvalue
	thermo["fO2low"], MeltInclusion_deep_thermo["logfO2"] = fO2_from_dQFM(MeltInclusion_deep_thermo)

	#Calculate all fugacity coefficient, logK, and K values and put them in the thermo dict
	CalcAllGammas(thermo)
	Calc_all_logK(thermo)
	logK_to_K(thermo)

	respeciate(thermo, 1) #takes args thermodynamic dict, new pressure in bars. Outputs all X's via Calc_all_Xs_from_fS method

	if thermo["solvefS2"] == 0:
		pass
	else:
		x_values_dict = {	"XCO": thermo["XCO"],
							"XCO2": thermo["XCO2"],
							"XH2": thermo["XH2"],
							"XH2O": thermo["XH2O"],
							"XH2S": thermo["XH2S"],
							"XO2": thermo["XO2"],
							"XS2": thermo["XS2"],
							"XSO2": thermo["XSO2"]}

		if all(value >= 0 for value in x_values_dict.values()):
			#Normalize X values
			Normalized_fluid = normalize_final_Xs(thermo)

			if verbose == True:
				#Print new X values at 1 bar
				print "\n"
				print "Sample: " + str(sample)
				print "\n"
				print "Mole fractions:"
				print("\n".join("{}: {}".format(k, v) for k, v in Normalized_fluid.items()))
				print "\n"
				print "Molar gas ratios:"
				print "SO2/CO2 = " + str(Normalized_fluid["XSO2"] / Normalized_fluid["XCO2"])
				print "H2S/SO2 = " + str(Normalized_fluid["XH2S"] / Normalized_fluid["XSO2"])
				print "\n"

				return Normalized_fluid
			else:
				return Normalized_fluid
		else:
			print "Some values less than 0"


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


###-------------MATCHING MODEL----------------###
H2Orange = np.arange(minH2O, maxH2O, H2Ostep)
CO2range = np.arange(minCO2, maxCO2, CO2step)
Srange = np.arange(minS, maxS, Sstep)
Prange = np.arange(minP, maxP, Pstep)

iternumber = 0
list_of_fluid_dicts = [] #create an empty list which we will append each calculated dict (containing normalized calc'd fluid comps) to
thermolist = [] #create an empty list which we will append each thermo dict to
wtlist = [] #create an empty list which we will append petrology data to for each run
for i in H2Orange:
	for j in CO2range:
		for k in Srange:
			#Make your thermo dict with the values in this iteration.....
			MeltInclusion_deep_thermo = 	{	"highPressure": 2000, #bars
												"temp": 1000, #celcius
												"logfO2": -10.804261, #absolute logfO2 value #POAS is QFM
												"X_Sum": 0 #don't change this value
											} #dict

			MeltInclusion_shallow_thermo = {	"press": 1,
												"temp": 950,
												"logfO2": -11.782762, #POAS is QFM
											}

			MeltComposition_deep_wt = 	{	"SiO2": 55.4, #POAS sample P80b from Cigolini et al (1991)
											"TiO2": 0.94,
											"Al2O3": 16.9,
											"FeOstar": 0, #EITHER put FeOstar or separate FeO and Fe2O3, NOT BOTH! 
											"Fe2O3": 5.26,
											"FeO": 3.49,
											"MnO": 0.16,
											"MgO": 3.21,
											"CaO": 7.17,
											"Na2O": 3.38,
											"K2O": 2.11,
											"P2O5": 0.26,
											"H2O": i,
											"CO2": j,
											"S": k
										}

			MeltComposition_shallow_wt = 	{	"SiO2": 55.4, #POAS sample P80b from Cigolini et al (1991)
												"TiO2": 0.94,
												"Al2O3": 16.9,
												"FeOstar": 0,
												"Fe2O3": 5.26,
												"FeO": 3.49,
												"MnO": 0.16,
												"MgO": 3.21,
												"CaO": 7.17,
												"Na2O": 3.38,
												"K2O": 2.11,
												"P2O5": 0.26,
												"H2O": 0,
												"CO2": 0,
												"S": 0
											}

			#Some precursors.....
			#Store temp in Kelvin
			MeltInclusion_deep_thermo["tempK"] = MeltInclusion_deep_thermo["temp"] + 273.15

			#CALCULATE FLUID FROM DIFFERENCED MELT INCLUSIONS#
			DegassedMI_deep_to_shallow = Difference_MI_volatiles(MeltComposition_deep_wt, MeltComposition_shallow_wt)

			#SET LOW PRESSURE
			MeltInclusion_deep_thermo["press"] = 1.0

			#CONVERT WT% TO MOL FRACTION
			DegassedMI_X = Convert_wt_to_molfrac(DegassedMI_deep_to_shallow)

			iternumber += 1
			print iternumber

			MeltInclusion_deep_thermo["XCO2tot"] = DegassedMI_X["XCO2fl"] #volatiles must be done separately and nornmalized to volatiles-only
			MeltInclusion_deep_thermo["XH2Otot"] = DegassedMI_X["XH2Ofl"]
			MeltInclusion_deep_thermo["XStot"] = DegassedMI_X["XSfl"]

			#run the model.....
			MeltInclusion_deep_thermo["dont_save_this_iteration"] = False
			Normd_fluid = TheModel(MeltInclusion_deep_thermo) #THIS IS THE THERMODYNAMIC MODEL RUNNING#
			if Normd_fluid is None:
				if verbose_errors == True:
					print "\n"
					print "Verbose errors:"
					print "Normd_fluid:"
					print Normd_fluid
					print "Thermo dict:"
					print MeltInclusion_deep_thermo
				else:
					pass
				raise ValueError, "Normd_fluid is None"
			else:
				pass
			if MeltInclusion_deep_thermo["dont_save_this_iteration"] == False:
				Normd_fluid["ratio_SO2_CO2"] = Normd_fluid["XSO2"] / Normd_fluid["XCO2"]
				Normd_fluid["ratio_H2S_SO2"] = Normd_fluid["XH2S"] / Normd_fluid["XSO2"]

				#save data externally to this loop
				list_of_fluid_dicts.append(Normd_fluid)
				MeltInclusion_deep_thermo["iternumber"] = iternumber
				thermolist.append(MeltInclusion_deep_thermo)
				wtlist.append(MeltComposition_deep_wt)

			else:
				pass

thermoFrame = pandas.DataFrame(thermolist) 
wtFrame = pandas.DataFrame(wtlist)
data = pandas.DataFrame(list_of_fluid_dicts)

#calculate all S as SO2 and all C as CO2
data["XSO2_star"] = data["XSO2"] + data["XH2S"]
data["XCO2_star"] = data["XCO2"] + data["XCO"] * (0.5) * 3
data["ratio_SO2star_CO2star"] = data["XSO2_star"] / data["XCO2_star"]

alldataFrame = pandas.concat([thermoFrame, wtFrame, data], axis=1)
#alldataFrame.reset_index(drop=True, inplace=True)

#Save this new data to an Excel spreadsheet
writer = pandas.ExcelWriter(filename + '_output.xlsx', engine='xlsxwriter') #Create a Pandas Excel writer using XlsxWriter as the engine.
data.to_excel(writer, sheet_name='Surface Gas Ratios')
thermoFrame.to_excel(writer, sheet_name='Thermo Data')
wtFrame.to_excel(writer, sheet_name='Petro Data')
alldataFrame.to_excel(writer, sheet_name='All Data')
writer.save() #Close the Pandas Excel writer and output the Excel file

##----------------------- PLOTTING ----------------------##
#Do some plotting, if the user called for it
list_of_Xs = ['XCO', 'XCO2', 'XH2', 'XH2O', 'XH2S', 'XO2', 'XS2', 'XSO2']
if save_plots == True:
	#Draw the figure
	fig, ax = plt.subplots(4,4)

	ax[0,1].set_title('Gas Ratio Visualizer')
	divider = make_axes_locatable(ax[0,0])
	
	#SO2/CO2 versus everything:
	XCO_SO2_CO2 = ax[0,0].plot(alldataFrame["XCO"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	XCO2_SO2_CO2 = ax[0,1].plot(alldataFrame["XCO2"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	XH2_SO2_CO2 = ax[0,2].plot(alldataFrame["XH2"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	XH2O_SO2_CO2 = ax[0,3].plot(alldataFrame["XH2O"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	XH2S_SO2_CO2 = ax[1,0].plot(alldataFrame["XH2S"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	XO2_SO2_CO2 = ax[1,1].plot(alldataFrame["XO2"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	XS2_SO2_CO2 = ax[1,2].plot(alldataFrame["XS2"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	XSO2_SO2_CO2 = ax[1,3].plot(alldataFrame["XSO2"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')

	#H2S/SO2 versus everything:
	XCO_H2S_SO2 = ax[2,0].plot(alldataFrame["XCO"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	XCO2_H2S_SO2 = ax[2,1].plot(alldataFrame["XCO2"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	XH2_H2S_SO2 = ax[2,2].plot(alldataFrame["XH2"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	XH2O_H2S_SO2 = ax[2,3].plot(alldataFrame["XH2O"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	XH2S_H2S_SO2 = ax[3,0].plot(alldataFrame["XH2S"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	XO2_H2S_SO2 = ax[3,1].plot(alldataFrame["XO2"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	XS2_H2S_SO2 = ax[3,2].plot(alldataFrame["XS2"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	XSO2_H2S_SO2 = ax[3,3].plot(alldataFrame["XSO2"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')

	#add horizontal lines to represent measured values
	SO2_CO2_min = eruption_SO2_CO2-eruption_SO2_CO2_range
	SO2_CO2_max = eruption_SO2_CO2+eruption_SO2_CO2_range
	H2S_SO2_min = eruption_H2S_SO2-eruption_H2S_SO2_range
	H2S_SO2_max = eruption_H2S_SO2+eruption_H2S_SO2_range

	#q values represent maximum and minimum values read in from GasData excel file (stored in GasData DataFrame)
	q_SO2_CO2_min = minvalgasdata1 
	q_SO2_CO2_max = maxvalgasdata1
	q_H2S_SO2_min = minvalgasdata2
	q_H2S_SO2_max = maxvalgasdata2

	for X in range(4):
		ax[0,X].axhspan(SO2_CO2_min, SO2_CO2_max, alpha=1, color='orange')
		ax[1,X].axhspan(SO2_CO2_min, SO2_CO2_max, alpha=1, color='orange')
		ax[2,X].axhspan(H2S_SO2_min, H2S_SO2_max, alpha=1, color='orange')
		ax[3,X].axhspan(H2S_SO2_min, H2S_SO2_max, alpha=1, color='orange')

		ax[0,X].axhspan(q_SO2_CO2_min, q_SO2_CO2_max, alpha=0.6, color='green')
		ax[1,X].axhspan(q_SO2_CO2_min, q_SO2_CO2_max, alpha=0.6, color='green')
		ax[2,X].axhspan(q_H2S_SO2_min, q_H2S_SO2_max, alpha=0.6, color='green')
		ax[3,X].axhspan(q_H2S_SO2_min, q_H2S_SO2_max, alpha=0.6, color='green')


	#Set labels for each subplot
	for i in range(4):
		ax[0,i].set_xlabel(list_of_Xs[i])
		ax[0,i].set_ylabel('SO2/CO2')

		ax[1,i].set_xlabel(list_of_Xs[i+4])
		ax[1,i].set_ylabel('SO2/CO2')

		ax[2,i].set_xlabel(list_of_Xs[i])
		ax[2,i].set_ylabel('H2S/SO2')

		ax[3,i].set_xlabel(list_of_Xs[i+4])
		ax[3,i].set_ylabel('H2S/SO2')

	#show the plot
	plt.show()

	#Second plot comparing calculated gas ratios to petrology (MI H2O, CO2, S wt percents) 
	fig2, ax2 = plt.subplots(2,3)

	ax2[0,2].set_title('Gas Ratios vs MI Petrology')
	divider2 = make_axes_locatable(ax2[0,0])

	#SO2/CO2 versus everything:
	wtH2O_SO2_CO2 = ax2[0,0].plot(alldataFrame["H2O"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	wtCO2_SO2_CO2 = ax2[0,1].plot(alldataFrame["CO2"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')
	wtS_SO2_CO2 = ax2[0,2].plot(alldataFrame["S"], alldataFrame["ratio_SO2_CO2"], color='red', marker='.', linestyle='none')

	#H2S/SO2 versus everything:
	wtH2O_H2S_SO2 = ax2[1,0].plot(alldataFrame["H2O"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	wtCO2_H2S_SO2 = ax2[1,1].plot(alldataFrame["CO2"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')
	wtS_H2S_SO2 = ax2[1,2].plot(alldataFrame["S"], alldataFrame["ratio_H2S_SO2"], color='blue', marker='.', linestyle='none')

	#add horizontal lines to represent measured values
	for X in range(3):
		ax2[0,X].axhspan(SO2_CO2_min, SO2_CO2_max, alpha=1, color='orange')
		ax2[1,X].axhspan(H2S_SO2_min, H2S_SO2_max, alpha=1, color='orange')

		ax2[0,X].axhspan(q_SO2_CO2_min, q_SO2_CO2_max, alpha=0.6, color='green')
		ax2[1,X].axhspan(q_H2S_SO2_min, q_H2S_SO2_max, alpha=0.6, color='green')

	#Set labels for each subplot
	list_of_vols = ["H2O", "CO2", "S"]
	for i in range(3):
		ax2[0,i].set_xlabel(list_of_vols[i])
		ax2[0,i].set_ylabel('SO2/CO2')

		ax2[1,i].set_xlabel(list_of_vols[i])
		ax2[1,i].set_ylabel('H2S/SO2')

	#show the plot
	plt.show()


	#make 100% stacked area chart to visualize mole fraction compositions
	#copy the data dataframe and drop some columns...
	data_area_plot = data.drop('X_Sum', axis=1) 
	data_area_plot = data_area_plot.drop('ratio_H2S_SO2', axis=1)
	data_area_plot = data_area_plot.drop('ratio_SO2_CO2', axis=1)
	data_area_plot = data_area_plot.drop('ratio_SO2star_CO2star', axis=1)
	data_area_plot = data_area_plot.drop('XSO2_star', axis=1)
	data_area_plot = data_area_plot.drop('XCO2_star', axis=1)


	# #sort by increasing CO2
	# data_area_plot = data_area_plot.sort_values(["XCO2", "XH2O"], ascending=True)
	
	#plot
	data_area_plot.plot.area()
	plt.show()

print '\n'
print 'SUCCESS! Saved file ' + filename + '_output.xlsx'

######-------Import timeseries gas data and match to calculated values-------#####
#Determines which dissolved volatile concentrations are required to reproduce the timeseries data.
positiveMatch = pandas.DataFrame() #Empty dataframe that stores values only when a match is positive for both ratios
alldataFormatch = pandas.DataFrame()
iternumber = 1
if MatchBoth == False:
	for index1, row1 in GasData.iterrows():
		print row1
		print iternumber
		iternumber += 1
		for index, row in alldataFrame.iterrows():
			matchH2S_SO2_min = row["ratio_H2S_SO2"] - tolerance * row["ratio_H2S_SO2"]
			matchH2S_SO2_max = row["ratio_H2S_SO2"] + tolerance * row["ratio_H2S_SO2"]
			matchSO2_CO2_min = row["ratio_SO2_CO2"] - tolerance * row["ratio_SO2_CO2"]
			matchSO2_CO2_max = row["ratio_SO2_CO2"] + tolerance * row["ratio_SO2_CO2"]

			if matchH2S_SO2_min <= row1["H2S/SO2"] <= matchH2S_SO2_max:
				if minvalgasdata1 <= row["ratio_SO2_CO2"] <= maxvalgasdata1:
					if minvalgasdata2 <= row["ratio_H2S_SO2"] <= maxvalgasdata2:
						positiveMatch = positiveMatch.append({	"Time": row1["Time"], 
																"H2S/SO2": row1["H2S/SO2"],
																"H2Omelt": row["H2O"],
																"CO2melt": row["CO2"],
																"Smelt": row["S"]}, ignore_index=True)
			alldataFormatch = alldataFormatch.append(row)

			if matchSO2_CO2_min <= row1["SO2/CO2"] <= matchSO2_CO2_max:
				if minvalgasdata1 <= row["ratio_SO2_CO2"] <= maxvalgasdata1:
					if minvalgasdata2 <= row["ratio_H2S_SO2"] <= maxvalgasdata2:
						positiveMatch = positiveMatch.append({	"Time": row1["Time"],
																"SO2/CO2": row1["SO2/CO2"],
																"H2Omelt": row["H2O"],
																"CO2melt": row["CO2"],
																"Smelt": row["S"]}, ignore_index=True)
						print row
						alldataFormatch = alldataFormatch.append(row)

if MatchBoth == True:
	for index1, row1 in GasData.iterrows():
		print row1
		print iternumber
		iternumber += 1
		for index, row in alldataFrame.iterrows():
			matchH2S_SO2_min = row["ratio_H2S_SO2"] - tolerance * row["ratio_H2S_SO2"]
			matchH2S_SO2_max = row["ratio_H2S_SO2"] + tolerance * row["ratio_H2S_SO2"]
			matchSO2_CO2_min = row["ratio_SO2_CO2"] - tolerance * row["ratio_SO2_CO2"]
			matchSO2_CO2_max = row["ratio_SO2_CO2"] + tolerance * row["ratio_SO2_CO2"]

			if matchH2S_SO2_min <= row1["H2S/SO2"] <= matchH2S_SO2_max and matchSO2_CO2_min <= row1["SO2/CO2"] <= matchSO2_CO2_max:
				if minvalgasdata1 <= row["ratio_SO2_CO2"] <= maxvalgasdata1 and minvalgasdata1 <= row["ratio_SO2_CO2"] <= maxvalgasdata1:
					if minvalgasdata2 <= row["ratio_H2S_SO2"] <= maxvalgasdata2 and minvalgasdata2 <= row["ratio_H2S_SO2"] <= maxvalgasdata2:
						positiveMatch = positiveMatch.append({	"Time": row1["Time"], 
																"H2S/SO2": row["ratio_H2S_SO2"],
																"SO2/CO2": row["ratio_SO2_CO2"],
																"H2Omelt": row["H2O"],
																"CO2melt": row["CO2"],
																"Smelt": row["S"]}, ignore_index=True)
						alldataFormatch = alldataFormatch.append({	"XH2Otot": row["XH2Otot"],
																	"XStot": row["XStot"],
																	"XCO2tot": row["XCO2tot"],
																	"press": row["press"],
																	"temp": row["temp"],
																	"tempK": row["tempK"]}, ignore_index=True)
print positiveMatch

#Save this new data to an Excel spreadsheet
writer2 = pandas.ExcelWriter(filename + '_Match.xlsx', engine='xlsxwriter') #Create a Pandas Excel writer using XlsxWriter as the engine.
positiveMatch.to_excel(writer2, sheet_name='All Matches')
GasData.to_excel(writer2, sheet_name='Gas Data')
alldataFormatch.to_excel(writer2, sheet_name='All Data')
writer2.save() #Close the Pandas Excel writer and output the Excel file

#plot user input gas data
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(GasData["Time"], GasData["SO2/CO2"], color='green', marker='.', linestyle='none')
axarr[0].plot(positiveMatch["Time"], positiveMatch["SO2/CO2"], color='yellow', marker='o', linestyle='none')
axarr[0].set_yscale('log')
axarr[0].set_title('Real Gas Data From Maarten')
axarr[1].plot(GasData["Time"], GasData["H2S/SO2"], color='blue', marker='.', linestyle='none')
axarr[1].plot(positiveMatch["Time"], positiveMatch["H2S/SO2"], color='yellow', marker='o', linestyle='none')
axarr[0].set_ylabel("SO2/CO2")
axarr[1].set_ylabel("H2S/SO2")
axarr[0].set_xlabel("Date")

#show the plots
plt.show()