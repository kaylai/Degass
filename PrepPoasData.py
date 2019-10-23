import pandas
import datetime
import matplotlib.pyplot as plt

#Get data and remove duplicate times (average if duplicate)
GasData = pandas.read_excel("MultiGasDataExcel_Fresh.xlsx")
GasData = GasData.groupby(["Time"], as_index=False).mean()
GasData = GasData.dropna()

#write data to new excel sheet
writer = pandas.ExcelWriter('PreppedPoasData.xlsx', engine='xlsxwriter')
GasData.to_excel(writer)
writer.save()

#plot the data
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(GasData["Time"], GasData["SO2/CO2"], color='green', marker='.', linestyle='none')
axarr[0].set_yscale('log')
axarr[0].set_title('Real Gas Data From Maarten')
axarr[1].plot(GasData["Time"], GasData["H2S/SO2"], color='blue', marker='.', linestyle='none')

#show the plot
plt.show()