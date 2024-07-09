#Written by Reuben Neate
#Processing WVR tipping curve data
#BEGIN-------------------------------------------------------------------------
#------------------------------------------------------------------------------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats

#------------------------------------------------------------------------------
#Retrieving the antenna temperature using eq 5 from the technical documentation
def RetrieveAntTemp(Pa, Ph, Pc, Th, Tc):
    #return ((Pa - Pc)*(Ph - Pc)/(Th - Tc)) + Tc
    g = (Th-Tc)/(Ph-Pc)
    o = Tc - (g*Pc)
    return(g*Pa) + o
#------------------------------------------------------------------------------
def DataExtract(data, Th_23, Th_31):
    angles = [0]*8
    temp_23 = [0]*8
    temp_31 = [0]*8
    for i in range(0,8):
         angles[i] = float(data[i][2])
         Tc = data[i][10]
         temp_23[i] = RetrieveAntTemp(data[i][8], data[i][7], data[i][6], Th_23, Tc)
         Tc = data[i][12]
         temp_31[i] = RetrieveAntTemp(data[i][5], data[i][4], data[i][3], Th_31, Tc)
    return angles, temp_23, temp_31        
#------------------------------------------------------------------------------
#Calculate airmass
def Airmass(A):
   return 1/np.sin(np.pi*A/180) 
#------------------------------------------------------------------------------
#Get the line characteristics
def LineFit(x_data, y_data):
    slope, intercept, r, p, std_err = stats.linregress(x_data, y_data)
    return slope, intercept
#------------------------------------------------------------------------------
#Returns the noise source hot temperature as from the tip curves
def TipCurveCalibrate(airmass, dataline, accuracy): 
    Th_23 = 700
    Th_31 = 700
    error_23 = 1000
    error_31 = 1000
    while((abs(error_23) > accuracy) or (abs(error_31) > accuracy)):    
        angles, temp_23, temp_31 = DataExtract(dataline, Th_23, Th_31)            
        if (abs(error_23) > accuracy): 
            # print(Th_23)
            slope_23, intercept_23 = LineFit(airmass, temp_23)
            error_23 = intercept_23 - 2.75
            Th_23 = Th_23 + error_23/10
        if (abs(error_31) > accuracy):
            slope_31, intercept_31 = LineFit(airmass, temp_31)
            error_31 = intercept_31 - 2.75
            Th_31 = Th_31 + error_31/10
    return Th_23, Th_31
#------------------------------------------------------------------------------
def PlotTipCurves(angles, temp_23, temp_31, title):
    airmass_2 = [0]*9
    airmass = [0]*8
    for i in range(8):
        airmass[i] = Airmass(angles[i])
        airmass_2[i+1] = Airmass(angles[i])
    slope_23, intercept_23 = LineFit(airmass, temp_23)
    slope_31, intercept_31 = LineFit(airmass, temp_31)
    line_23 = [0]*9
    line_31 = [0]*9
    for i in range(len(airmass_2)):
        line_23[i] = (airmass_2[i] * slope_23) + intercept_23
        line_31[i] = (airmass_2[i] * slope_31) + intercept_31
    print("23 GHz 0 airmass intercept of ", str(intercept_23), " K")
    print("31 GHz 0 airmass intercept of ", str(intercept_31), " K")
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle(title)
    ax1.plot(airmass_2, line_23)
    ax1.plot(airmass, temp_23, 'x')
    ax1.set_title("23 GHz")
    ax1.set_xlabel("Airmass")
    ax1.set_ylabel("Brightness temperature (K)")
    ax1.grid()
    ax2.plot(airmass_2, line_31)
    ax2.plot(airmass, temp_31, 'x')
    ax2.set_title("31 Ghz")
    ax2.set_xlabel("Airmass")
    ax2.grid()
    plt.show()
    
    return
#------------------------------------------------------------------------------
def ProcessWVR3(file):
    """
    A function that process a datafile from WVR3
    
    Args:
        file (string): Name of data file to be processed
        
    Returns:
        date (string): date the data was taken on
        time (list of string): time the tip curve was started
        angles (list of list of float): Elevation angle in degrees
        temp_23 (list of list of float): Calibrated 23.8 GHz brightness temperature corresponding to the given angles
        temp_31 (list of list of float): Calibrated 31.6 GHz brightness temperature corresponding to the given angles
    """
    #Reading in the data
    with open(file) as data_file:
        data = []
        i = 0
        for line in data_file:
            data_1 = line.split(",")
            data_2 = [0]*len(data_1)
            for i in range(len(data_1)):
                if (i == 1) or (i == 2):
                    data_2[i] = data_1[i]
                else:
                    data_2[i] = (float(data_1[i]))
            data.append(data_2)   
#------------------------------------------------------------------------------
    #Check that the data starts at the an angle above 85 degrees, close enough to Zenith
    while (float(data[0][2]) < 85):
        del(data[0])
        
    l = len(data)//8
    p = int(8*((len(data)/8) - l))
    
    t23 = [] #list to store the 23.8 GHz brightness temeprature output
    t31 = [] #list to store the 31.6 GHz brightness temeprature output
    a = [] #list to store the angles measured at
    date = data[0][1].split(" ")[1]
    time = []
    for n in range(0,l,8):
        time.append(data[n][1].split(" ")[2])
    
    for n in range(l): #Loop through a and process all of the data
#------------------------------------------------------------------------------
        #Reading the first 8 lines because that is a full tip curve
        s = (8*n)
        Th_nom = 700 # nominal or starting guess for hot load temperature
        angles, temp_23, temp_31 = DataExtract(data[s:s+8], Th_nom, Th_nom)
        
        #Creating an array of arrmass form the angles
        airmass = [0]*8
        for i in range(8):
            airmass[i] = Airmass(angles[i])
            
        # if n == 10: #For testing purposes you can plot the tip curve data
            # print("Raw tip curve data")
            # print("Expected hot load temperature: 700 K")
            # #PlotTipCurves(angles, temp_23, temp_31, "Uncorrected Tip curves")
            # print("")
#------------------------------------------------------------------------------
        #Apply tipping curve calibration to correct the data
        accuracy = 0.001  
        Th_23, Th_31 = TipCurveCalibrate(airmass, data[s:s+8], accuracy)
        angles, temp_23, temp_31 = DataExtract(data[s:s+8], Th_23, Th_31) 
    
        # if n == 10: #For testing purposes you can plot the tip curve data
        #     print("Calibrated tip curve data")  
        #     print("23 GHz hot load temperature: ", Th_23, " K")
        #     print("31 GHz hot load temperature: ", Th_31, " K")
        #     PlotTipCurves(angles, temp_23, temp_31, "Corrected Tip curves")
        #     print("")
        
        a.append(angles)
        t23.append(temp_23)
        t31.append(temp_31)
    
    return date, time, a, t23, t31
            
         
#------------------------------------------------------------------------------







