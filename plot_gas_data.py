import numpy as np
import pandas
import matplotlib.pyplot as plt 

###For algorithm to match with surface gas data
gasdatafilename = "MultiGasDataExcel.xlsx" #filename of excel file with surface gas data
tolerance = 0.20 #percent allowable tolerance around calcaulted values to find a match with observed data. 0.05 = 5%

GasData = pandas.read_excel(gasdatafilename) #Read in gas data

#Sort GasData by data type to match Maarten's plot...
LakeGasFixed = GasData[GasData['ID_Loc'].str.contains('Lake')]
DomeMobileMG = GasData[GasData["ID_Loc"].str.contains('Dome')]
Drone = GasData[GasData["ID_Loc"].str.contains('Drone')]
OVSI_MG1 = GasData[GasData["ID_Loc"].str.contains('Rim')]
MainVent = GasData[GasData["ID_Loc"].str.contains('MainVent')]
DryLakeVents = GasData[GasData["ID_Loc"].str.contains('unknown')]

print DomeMobileMG

#plot user input gas data
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(GasData["Time"], GasData["SO2/CO2"], color='blue', marker='.', linestyle='none')
axarr[0].plot(DomeMobileMG["Time"], DomeMobileMG["SO2/CO2"], mfc='red', mec='black', marker='^', linestyle='none', markersize=6)
axarr[0].plot(LakeGasFixed["Time"], LakeGasFixed["SO2/CO2"], mfc='yellow', mec='black', marker='o', linestyle='none')
axarr[0].plot(Drone["Time"], Drone["SO2/CO2"], color='orange', mec='black', marker='D', linestyle='none', markersize=6)
axarr[0].plot(MainVent["Time"], MainVent["SO2/CO2"], color='red', mec='black', marker='^', linestyle='none', markersize=6)
axarr[0].plot(OVSI_MG1["Time"], OVSI_MG1["SO2/CO2"], color='#00FF00', mec='black', marker='o', linestyle='none')
axarr[0].set_yscale('log')
axarr[0].set_title('Real Gas Data From Maarten')
axarr[1].plot(GasData["Time"], GasData["H2S/SO2"], color='blue', marker='.', linestyle='none')
axarr[0].set_ylabel("SO2/CO2")
axarr[1].set_ylabel("H2S/SO2")
axarr[0].set_xlabel("Date")

#show the plot
plt.show()

#plot gas ratios against each other
f, ax = plt.subplots()
ax.plot(GasData["SO2/CO2"], GasData["H2S/SO2"], color='purple', marker='.', linestyle='none')
ax.set_title('Gas Ratios')
plt.xlabel('SO2/CO2')
plt.ylabel('H2S/SO2')

#show the plot
plt.show()

#plot user input gas data
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(GasData["Time"], GasData["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].plot(DomeMobileMG["Time"], DomeMobileMG["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].plot(LakeGasFixed["Time"], LakeGasFixed["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].plot(Drone["Time"], Drone["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].plot(MainVent["Time"], MainVent["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].plot(OVSI_MG1["Time"], OVSI_MG1["SO2/CO2"], color='black', marker='.', linestyle='none')
axarr[0].set_yscale('log')
axarr[0].set_title('Real Gas Data From Maarten')
axarr[1].plot(GasData["Time"], GasData["H2S/SO2"], color='black', marker='.', linestyle='none')
axarr[0].set_ylabel("SO2/CO2")
axarr[1].set_ylabel("H2S/SO2")
axarr[0].set_xlabel("Date")

#show the plot
plt.show()