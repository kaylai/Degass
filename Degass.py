from math import *

##-----------------USER INPUTS---------------##
#MELT INCLUSION WT% AS (H2O, CO2, S)
#Temp in degrees celsius
#Press in bars

sample = "erebus"

if sample == "poas":
	MeltInclusion_deep_thermo = 	{	"press": 2000, #bars
										"temp": 1000, #celcius
										"logfO2": -10.804261, #absolute logfO2 value #POAS is QFM
										"XCO2": 0.1701, #calculated with magmasat
										"X_Sum": 0 #don't change this value
									} #dict

	MeltInclusion_shallow_thermo = {	"press": 1,
									"temp": 950,
									"logfO2": -11.782762, #POAS is QFM
									"XCO2": 0.016 #calculated with magmasat
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
								"H2O": 5.0,
								"CO2": 0.2,
								"S": 0.2
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
										"H2O": 0.4,
										"CO2": 0.01,
										"S": 0.01 
									}

if sample == "erebus":
	MeltInclusion_deep_thermo = 	{	"press": 4445, #bars
										"temp": 1100, #celcius
										"logfO2": -7.63, 
										"XCO2": 0.93, 
										"X_Sum": 0 #don't change this value
									} #dict
	
	MeltInclusion_shallow_thermo = {	"press": 3582,
									"temp": 1081,
									"logfO2": -9.99
								}

	MeltComposition_deep_wt = 	{	"SiO2": 41.70, #POAS sample P80b from Cigolini et al (1991)
								"TiO2": 4.21,
								"Al2O3": 14.80,
								"FeOstar": 0, #EITHER put FeOstar or separate FeO and Fe2O3, NOT BOTH! 
								"Fe2O3": 5.72,
								"FeO": 8.05,
								"MnO": 0.16,
								"MgO": 5.95,
								"CaO": 13.07,
								"Na2O": 3.74,
								"K2O": 1.65,
								"P2O5": 0.94,
								"H2O": 1.50,
								"CO2": 0.5537,
								"S": 0.2166
							}

	MeltComposition_shallow_wt = 	{	"SiO2": 51.67, #POAS sample P80b from Cigolini et al (1991)
										"TiO2": 1.84,
										"Al2O3": 19.49,
										"FeOstar": 0,
										"Fe2O3": 1.56,
										"FeO": 6.64,
										"MnO": 0.21,
										"MgO": 1.33,
										"CaO": 4.02,
										"Na2O": 7.73,
										"K2O": 4.96,
										"P2O5": 0.55,
										"H2O": 0.28,
										"CO2": 0.1638,
										"S": 0.0918 
									}


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

##----------------THE MODEL-----------------##
#Store temp in Kelvin
MeltInclusion_deep_thermo["tempK"] = MeltInclusion_deep_thermo["temp"] + 273.15

#CALCULATE FLUID FROM DIFFERENCED MELT INCLUSIONS#
DegassedMI_deep_to_shallow = Difference_MI_volatiles(MeltComposition_deep_wt, MeltComposition_shallow_wt)

#CONVERT WT% TO MOL FRACTION
DegassedMI_X = Convert_wt_to_molfrac(DegassedMI_deep_to_shallow)
MeltInclusion_deep_thermo["XCO2"] = DegassedMI_X["XCO2fl"] #volatiles must be done separately and nornmalized to volatiles-only
MeltInclusion_deep_thermo["XH2O"] = DegassedMI_X["XH2Ofl"]
MeltInclusion_deep_thermo["XStot"] = DegassedMI_X["XSfl"]

#CALCULATE FO2 FROM LOGFO2
MeltInclusion_deep_thermo["fO2"] = 10**(MeltInclusion_deep_thermo["logfO2"])

##Speciate Degassed MI##
#Calculate fugacity coefficients and put them in the thermo dict
fugacity_coefficients_list = ['gammaCO', 'gammaCO2', 'gammaH2', 'gammaH2O', 'gammaH2S', 'gammaO2', 'gammaS2', 'gammaSO2'] #List of all fugacity coefficients
CP_list = [CPCO, CPCO2, CPH2, CPH2O, CPH2S, CPO2, CPS2, CPSO2] #List of all critical parameter dicts

CalcAllGammas(MeltInclusion_deep_thermo)

#Calculate all logK and K values and put them in the thermo dict
Calc_all_logK(MeltInclusion_deep_thermo)
logK_to_K(MeltInclusion_deep_thermo)

#Calculate all pure fugacities
Calc_All_Pure_fs(MeltInclusion_deep_thermo)

##CO2 FIRST CALC##
#Calculate fCO2
MeltInclusion_deep_thermo["fCO2"] = Calc_fCO2_from_XCO2(MeltInclusion_deep_thermo)


##O2##
#Calculate XO2
MeltInclusion_deep_thermo["XO2"] = Calc_XO2(MeltInclusion_deep_thermo)


##H2O##
#Calculate fH2O at high P using Moore et al 1998
#MeltInclusion_deep_thermo["fH2O"] = Calc_fH2O_Moore1998(DegassedMI_X, MeltInclusion_deep_thermo)
#Calculate fH2O from XH2O
MeltInclusion_deep_thermo["fH2O"] = Calc_fH2O_from_XH2O(MeltInclusion_deep_thermo)

#Calculate XH2O
#MeltInclusion_deep_thermo["XH2O"] = Calc_XH2O(MeltInclusion_deep_thermo)


##H2##
#Calculate fH2 at high P based on given fO2 and calculated fH2O
MeltInclusion_deep_thermo["fH2"] = Calc_fH2(MeltInclusion_deep_thermo)

#Calculate XH2
MeltInclusion_deep_thermo["XH2"] = Calc_XH2(MeltInclusion_deep_thermo)

#while MeltInclusion_deep_thermo["X_Sum"] < 0.999 or MeltInclusion_deep_thermo["X_Sum"] > 1.001:
##CO##
#Calcualte fCO
MeltInclusion_deep_thermo["fCO"] = Calc_fCO(MeltInclusion_deep_thermo)

#Calculate XCO
MeltInclusion_deep_thermo["XCO"] = Calc_XCO(MeltInclusion_deep_thermo)


##CO2 SECOND CALC##
#Check XCO2 is not too large
Check_XCO2_value(MeltInclusion_deep_thermo)


##SULFUR##
#Calculate non-S Partial Pressures
MeltInclusion_deep_thermo["PCO2"] = Calc_PCO2(MeltInclusion_deep_thermo)
MeltInclusion_deep_thermo["PH2O"] = Calc_PH2O(MeltInclusion_deep_thermo)
MeltInclusion_deep_thermo["PH2"] = Calc_PH2(MeltInclusion_deep_thermo)
MeltInclusion_deep_thermo["PCO"] = Calc_PCO(MeltInclusion_deep_thermo)
MeltInclusion_deep_thermo["PO2"] = Calc_PO2(MeltInclusion_deep_thermo)

#Calculate PStot
#MeltInclusion_deep_thermo["PStot"] = Calc_PStot(MeltInclusion_deep_thermo)
MeltInclusion_deep_thermo["PStot"] = MeltInclusion_deep_thermo["press"] * MeltInclusion_deep_thermo["XStot"]

#Calculate fS2 (using equation 7 from Iacovino (2015) EPSL)
MeltInclusion_deep_thermo["fS2"] = Calc_fS2(MeltInclusion_deep_thermo)

#Calculate fSO2
MeltInclusion_deep_thermo["fSO2"] = Calc_fSO2(MeltInclusion_deep_thermo)

#Calculate fH2S
MeltInclusion_deep_thermo["fH2S"] = Calc_fH2S(MeltInclusion_deep_thermo)

#Calculate XS2
MeltInclusion_deep_thermo["XS2"] = Calc_XS2(MeltInclusion_deep_thermo)


#Calculate Partial Pressures
Calc_all_PartialPressures(MeltInclusion_deep_thermo)

#Calculate fluid ratios
Calc_all_fluid_ratios(MeltInclusion_deep_thermo, DegassedMI_X)

#Normalize final mole fraction fluid values
Normalized_fluid = normalize_final_Xs(MeltInclusion_deep_thermo)

print_all_Xs(MeltInclusion_deep_thermo)
print Normalized_fluid

print MeltInclusion_deep_thermo["X_Sum"]


print("{" + "\n".join("{}: {}".format(k, v) for k, v in MeltInclusion_deep_thermo.items()) + "}")
#print("\n")
#print("current problem: ")


#dec15,2017 - maybe need totally different approach to speciate fluid, which uses XH2O, XCO2, and XStot inputs as total H, total C, total S, then speciates each of those.



