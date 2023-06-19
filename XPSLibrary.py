#XPS Plots

#Libraries needed
import numpy as np
from scipy import stats
import pandas as pd


#Given the path to the file as folde and file name, it reutrns the temperature ramp as a Pandas dataframe.
#The time is in a strange format and needs to be modified depending on the format chosen during the measurments
def ReadTemperatureData(FolderTR,FileTR):
    TemperatureRamp=pd.read_csv(FolderTR+"/"+FileTR+".csv") #Read the temperature ramp data
    #give meaningful names to the dataset columns
    TemperatureRamp.rename(columns={"Date/Time": "Time", "COM9.ID001-3216.INPUT.PVInValue": "Temperature","COM9.ID001-3216.SP.SP1": "Setpoint","COM9.ID001-3216.INPUT.PVDot": "Unknown"}, inplace=True)
    return TemperatureRamp

def AdjustTime(TemperatureRamp, unit):
    #Time is in a strange format. Subtracting the first time you get the time starting from zero.
    #Then, I convert the time in hours by multplying for the unit time I select while acquiring the data that was 24h
    TimeUnit=unit #hours
    TemperatureRamp['Time']=(TemperatureRamp['Time']-TemperatureRamp['Time'][0])*TimeUnit
    return TemperatureRamp
#counts the lines in a file
def Count_Lines(File):
    file = open(File, "r")
    lines = file.read().splitlines()
    return len(lines)

#counts the number of spectra in an XPS file by counitng how many times the line "Names\t\t" repeats itself
def Count_DataPerSpectrum(File):
    file = open(File, "r")
    lines = file.read().splitlines()
    NumberOfSpectra = 0
    for line in lines:
        print(line[0:4])
        if line[0:4] == ("Name"):
            NumberOfSpectra += 1
    file.close()
    DataPerSpectrum=len(lines)//NumberOfSpectra
    return DataPerSpectrum

#reads the XPS Data and gives meaninigful name to the columns. The colums must be 4
def ReadXPSData(FolderTR, FileXPS):
    HeaderColumns=8
    DataPerSpectrum=Count_DataPerSpectrum(FolderTR+"/"+FileXPS)
    XPS=pd.read_table(FolderTR+"/"+FileXPS, skiprows=lambda x: x%DataPerSpectrum<HeaderColumns, header=None)
    XPS=XPS.drop(columns=XPS.columns[2])
    XPS.columns=["Kinetic energy", "Counts", "Binding Energy", "CPS"]
    return XPS

#reads the XPS Data with the fitting functions, the name of the columns must the passed to the function
def ReadXPSDataFitted(FolderTR, FileXPS, HeaderColumns, ColumnsNames):
    DataPerSpectrum=Count_DataPerSpectrum(FolderTR+"/"+FileXPS)
    XPS=pd.read_table(FolderTR+"/"+FileXPS, skiprows=lambda x: x%DataPerSpectrum<HeaderColumns, header=None, sep="\t")
    XPS=XPS.dropna(axis='columns')
    XPS.columns=ColumnsNames
    return XPS

#Add the temperature to the temperature column at which the acquisition started
def SetTemperatureColumn(XPS, NumberOfSpectra, regression):
    for i in range(0, NumberOfSpectra):
        XPS["Initial Temperature"][i*NumericDataPerSpectrum:(i+1)*NumericDataPerSpectrum]=regression.intercept+regression.slope
    return XPS

#Returns a vector with the temperature at which each spectrum was acquired
def TemperatureVector(regression, time, NumberOfSpectra):
    t=np.arange(NumberOfSpectra+1)
    Duration=(np.max(time)-np.min(time))/NumberOfSpectra
    Temperature=regression.intercept+regression.slope*Duration*t
    return Temperature
    
#It computer the mean of a region of CPS (greater than a setpoint named lower energy or lower than a setpoint named upper energy) 
#and subtract it to to the data for each spectrum so that all the spectra cross the zero line.
def Normalize(XPS, NumberOfSpectra, NumberOfData):  
    
    x_len=XPS.shape[0]//NumberOfSpectra
    #reshape the data to compute the mean for each spectrum
    Means=np.mean(XPS["CPS"].values.reshape(NumberOfSpectra,x_len)[0:NumberOfData],axis=0)
    #subtratct the mean to each spectrum 
    for i in range (0, NumberOfSpectra):
        XPS["CPS"].values.reshape(NumberOfSpectra,x_len)[i]-= Means[i]
    
    return XPS
    
#For each spectrum, it computes the linear background and subtractit to the data
#To do so, it uses two doublets (x,y),
#they are the mean of the Point number of points at the right of the left edge and at the left of the rigth edge
def RemoveLinearBackgorund(XPS, Points, NumberOfSpectra):
    x_len=XPS.shape[0]//NumberOfSpectra

    for i in range (0, NumberOfSpectra):
        #for each spectrum compute the doublets (x,y) that will be used to calculate the slope m and intercept q
        y1=np.mean(XPS["CPS"].values.reshape(NumberOfSpectra,x_len)[i][0:Points])
        y2=np.mean(XPS["CPS"].values.reshape(NumberOfSpectra,x_len)[i][x_len-Points:x_len])

        x1=np.mean(XPS["Binding Energy"].values.reshape(NumberOfSpectra,x_len)[i][0:Points])
        x2=np.mean(XPS["Binding Energy"].values.reshape(NumberOfSpectra,x_len)[i][x_len-Points:x_len])

        #compute slope and intercept
        m=(y2-y1)/(x2-x1)
        q=y2-m*x2

        #subtract to the CPS the as computed linear background
        XPS["CPS"].values.reshape(NumberOfSpectra,x_len)[i] -= m*XPS["Binding Energy"].values.reshape(NumberOfSpectra,x_len)[i] + q
        
    return XPS

def NormalizeAU(Counts,RefMax):
    Counts= Counts/RefMax #rescale the counts
    return Counts

##########################################################################################################################################
#Report display and analysis
import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact, interactive, fixed, interact_manual
from IPython.display import display
import pandas as pd

#Read text file and divide in into a an array whose entries are the columns
def FileToArray(filename):
    file = open(filename + ".txt", "r")
    lines = file.read().splitlines() #puts the file into an array
    #lines
    return lines

#Since there are some missing data in the file, it is not possible to open the file as a astropy table. So the lines must be filled
def FillArray(lines):
    text= ''
    for i in range(1, len(lines)):

        divided_line=lines[i].split("\t")
        #print('Line', divided_line)

        if divided_line[0] != '':
            text= divided_line[0]
            #print('Text to write:', text)

        if divided_line[0] == '':
            divided_line[0] = text
            #print('New line:', divided_line)
        #print('\n')
        lines[i]=divided_line
    return lines

#Printing the correct lines in a new file
def PrintFilledFile(filename, lines):
    
    file = open(filename+"_filled"+".txt", "w")

    #print title
    for w in lines[0]:
        file.write(str(w))
    file.write("\n")

    #print data
    for i in range(1,len(lines)):
        l=lines[i]
        for w in l:
            file.write(str(w))
            if w != l[len(l)-1]:
                file.write("\t")
        file.write("\n")
    file.close()