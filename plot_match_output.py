import numpy as np
import pandas
import matplotlib.pyplot as plt 
import matplotlib
import datetime

###Bring in file to plot output from Matching alorithm in Degass.py
matchdatafilename = "Poas_PrimitiveBasalt_model_Match.xlsx" #filename of excel file with surface gas data
gasdatafilename = "MultiGasDataExcel_Fresh.xlsx" #filename of excel file with surface gas data

MatchData = pandas.read_excel(matchdatafilename) #Read in match data
GasData = pandas.read_excel(gasdatafilename) #Read in gas data

#Make a new datafame only containing rows where both SO2/CO2 and H2S/SO2 matched
#TODO

#Only take rows where SO2/CO2 is finite...
MatchDataSO2CO2 = MatchData[pandas.notnull(MatchData["SO2/CO2"])]

#Only take rows where H2S/SO2 is finite...
MatchDataH2SSO2 = MatchData[pandas.notnull(MatchData["H2S/SO2"])]

########################### PRINT SOME INTERESTING STUFF ###########################
print "Unique Ratios:"
print MatchData.groupby(["SO2/CO2", "H2S/SO2"]).count()
print "Minimum SO2/CO2 = " + str(min(MatchData["SO2/CO2"]))
print "Maximum SO2/CO2 = " + str(max(MatchData["SO2/CO2"]))
print "Minimum H2S/SO2 = " + str(min(MatchData["H2S/SO2"]))
print "Maximum H2S/SO2 = " + str(max(MatchData["H2S/SO2"]))

############################ PLOTTING ##############################

####plot for paper...

#format tick marks to show MM/YY
locator = matplotlib.dates.MonthLocator((1, 7))
formatter = matplotlib.dates.DateFormatter('%m/%y')

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)

#plot real gas data...
axarr[0].plot(GasData["Time"], GasData["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[1].plot(GasData["Time"], GasData["H2S/SO2"], color='black', marker='.', linestyle='none')


#overlay match data...
axarr[0].plot(MatchData["Time"], MatchData["SO2/CO2"], mfc='yellow', mec='black', marker='o', linestyle='none', markersize=6)
axarr[1].plot(MatchData["Time"], MatchData["H2S/SO2"], mfc='yellow', mec='black', marker='o', linestyle='none', markersize=6)

#overlay bars corresponding to "purely magmatic" gas ratios
axarr[0].axhspan(min(MatchData["SO2/CO2"]), max(MatchData["SO2/CO2"]), alpha=0.7, color='green')
axarr[1].axhspan(min(MatchData["H2S/SO2"]), max(MatchData["H2S/SO2"]), alpha=0.7, color='green')

#overlay vertical bars corresponding to phases in Maarten's paper...
axarr[0].axvspan(datetime.datetime(2012, 12, 1), datetime.datetime(2014, 6, 1), alpha=0.2, color='blue') #phase 1
axarr[0].axvspan(datetime.datetime(2016, 6, 1), datetime.datetime(2016, 9, 1), alpha=0.2, color='blue') #phase 2
axarr[0].axvspan(datetime.datetime(2017, 4, 1), datetime.datetime(2017, 9, 1), alpha=0.2, color='blue') #phase 3

axarr[1].axvspan(datetime.datetime(2012, 12, 1), datetime.datetime(2014, 6, 1), alpha=0.2, color='blue') #phase 1
axarr[1].axvspan(datetime.datetime(2016, 6, 1), datetime.datetime(2016, 9, 1), alpha=0.2, color='blue') #phase 2
axarr[1].axvspan(datetime.datetime(2017, 4, 1), datetime.datetime(2017, 9, 1), alpha=0.2, color='blue') #phase 3

#set labels, styles, ticks, etc
axarr[0].set_yscale('log')
axarr[0].set_ylabel("SO2/CO2")
axarr[1].set_ylabel("H2S/SO2")
axarr[0].set_xlim(datetime.datetime(2013,1,1), datetime.datetime(2018,1,1))
axarr[1].set_ylim(0,5)
axarr[0].xaxis.set(major_formatter=formatter, major_locator=locator)
axarr[1].xaxis.set(major_formatter=formatter, major_locator=locator)
f.autofmt_xdate()

#show the plot
plt.show()

########################## PLOTS THAT ARE TURNED OFF #####################
#Uncomment plt.show() to plot these...


####plot gas data with Matches overlain and thermo data over time for all matches
#Sort GasData by data type to match Maarten's plot...
LakeGasFixed = GasData[GasData['ID_Loc'].str.contains('Lake')]
DomeMobileMG = GasData[GasData["ID_Loc"].str.contains('Dome')]
Drone = GasData[GasData["ID_Loc"].str.contains('Drone')]
OVSI_MG1 = GasData[GasData["ID_Loc"].str.contains('Rim')]
MainVent = GasData[GasData["ID_Loc"].str.contains('MainVent')]
DryLakeVents = GasData[GasData["ID_Loc"].str.contains('unknown')]

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(5, sharex=True)

#plot real gas data...
axarr[0].plot(GasData["Time"], GasData["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].set_yscale('log')
axarr[0].set_title('Real Gas Data From Maarten')
axarr[1].plot(GasData["Time"], GasData["H2S/SO2"], color='black', marker='.', linestyle='none')
axarr[0].set_ylabel("SO2/CO2")
axarr[1].set_ylabel("H2S/SO2")
axarr[0].set_xlabel("Date")

#overlay match data...
axarr[0].plot(MatchData["Time"], MatchData["SO2/CO2"], mfc='yellow', mec='black', marker='o', linestyle='none', markersize=6)
axarr[1].plot(MatchData["Time"], MatchData["H2S/SO2"], mfc='yellow', mec='black', marker='o', linestyle='none', markersize=6)

#plot thermo data over time for matches...
axarr[2].plot(MatchData["Time"], MatchData["H2Omelt"], color='blue', marker='o', linestyle='none', markersize=6)
axarr[3].plot(MatchData["Time"], MatchData["CO2melt"], color='red', marker='o', linestyle='none', markersize=6)
axarr[4].plot(MatchData["Time"], MatchData["Smelt"], color='green', marker='o', linestyle='none', markersize=6)

#overlay bars corresponding to "purely magmatic" gas ratios
axarr[0].axhspan(min(MatchData["SO2/CO2"]), max(MatchData["SO2/CO2"]), alpha=0.9, color='gray')
axarr[1].axhspan(min(MatchData["H2S/SO2"]), max(MatchData["H2S/SO2"]), alpha=0.9, color='gray')


#show the plot
plt.show()

####plot gas data with Matches overlain and thermo data over time for only SO2/CO2 matches
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(4, sharex=True)

#plot real gas data...
axarr[0].plot(GasData["Time"], GasData["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].set_yscale('log')

#overlay match data...
axarr[0].plot(MatchDataSO2CO2["Time"], MatchDataSO2CO2["SO2/CO2"], mfc='yellow', mec='black', marker='o', linestyle='none', markersize=6)

#plot thermo data over time for matches...
axarr[1].plot(MatchDataSO2CO2["Time"], MatchDataSO2CO2["H2Omelt"], color='blue', marker='o', linestyle='none', markersize=6)
axarr[2].plot(MatchDataSO2CO2["Time"], MatchDataSO2CO2["CO2melt"], color='red', marker='o', linestyle='none', markersize=6)
axarr[3].plot(MatchDataSO2CO2["Time"], MatchDataSO2CO2["Smelt"], color='green', marker='o', linestyle='none', markersize=6)

#set labels
axarr[0].set_title('Real Gas Data From Maarten')
axarr[0].set_ylabel("SO2/CO2")
axarr[0].set_xlabel("Date")
axarr[1].set_ylabel("H2O melt")
axarr[2].set_ylabel("CO2 melt")
axarr[3].set_ylabel("S melt")

#show the plot
#plt.show()

####plot gas data with Matches overlain and thermo data over time for only H2S/SO2 matches
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(4, sharex=True)

#plot real gas data...
axarr[0].plot(GasData["Time"], GasData["H2S/SO2"], color='black', marker='.', linestyle='none')

#overlay match data...
axarr[0].plot(MatchDataH2SSO2["Time"], MatchDataH2SSO2["H2S/SO2"], mfc='yellow', mec='black', marker='o', linestyle='none', markersize=6)

#plot thermo data over time for matches...
axarr[1].plot(MatchDataH2SSO2["Time"], MatchDataH2SSO2["H2Omelt"], color='blue', marker='o', linestyle='none', markersize=6)
axarr[2].plot(MatchDataH2SSO2["Time"], MatchDataH2SSO2["CO2melt"], color='red', marker='o', linestyle='none', markersize=6)
axarr[3].plot(MatchDataH2SSO2["Time"], MatchDataH2SSO2["Smelt"], color='green', marker='o', linestyle='none', markersize=6)

#set labels
axarr[0].set_title('Real Gas Data From Maarten')
axarr[0].set_ylabel("H2S/SO2")
axarr[0].set_xlabel("Date")
axarr[1].set_ylabel("H2O melt")
axarr[2].set_ylabel("CO2 melt")
axarr[3].set_ylabel("S melt")

#show the plot
#plt.show()

# f, axarr = plt.subplots(2)

# axarr[0].plot(MatchData["Smelt"], MatchData["SO2/CO2"], color='green', marker='o', linestyle='none',  markersize=6)
# axarr[1].plot(MatchData["CO2melt"], MatchData["SO2/CO2"], color='red', marker='o', linestyle='none', markersize=6)

# plt.show()
