import numpy as np
import pandas
import matplotlib.pyplot as plt 

###Bring in file to plot output from Matching alorithm in Degass.py
matchdatafilename = "Poas_model_Match.xlsx" #filename of excel file with surface gas data
gasdatafilename = "MultiGasDataExcel.xlsx" #filename of excel file with surface gas data

MatchData = pandas.read_excel(matchdatafilename) #Read in match data
GasData = pandas.read_excel(gasdatafilename) #Read in gas data