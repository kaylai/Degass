import pandas

data = pandas.read_excel("PoasGasData.xlsx")

datetimevalue = 0
for index, row in data.iterrows():
	currentdatetimevalue = data["Time"]
	if currentdatetimevalue == datetimevalue:
		


	datetimevalue = data["Time"]