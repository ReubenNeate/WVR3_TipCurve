import pandas as pd
from glob import glob
import numpy as np
from scipy import stats
from sys import argv
import argparse
#import matplotlib.pyplot as plt


""" Requires path to file or to files; An argument with file path can be passed when running file

* Drop all radiometer csv files of each month in one folder/directory and provide path to analyze all data

 Example:


python TipCurves_dataframes.py -p '/Users/thecuriousvambo/22_rad/all_unprocces_22/*.csv'

or for one file:

python TipCurves_dataframes.py -p 'T240312.csv'

where -p is the arg to path file. Note the strings qoutes around the path of file or files

or an absolute path can be passed to line 346, and the one can then run the program without -p argument as:

python TipCurves_dataframes.py


"""


## this function is for secant of angle

def sec(theta):
    res = 1/np.cos(theta)
    return res


## this function is to calculate opacity using three brightness temp


def opacity_mean_temp(zenith_1, zenith_2, zenith_3):


    chan = [23,31]
    opac = []

    mean_temp = []

    for f in chan:
        ## calculate the opacity &

        a1 = list(zenith_1['zenith_angle'])[0]
        a2 = list(zenith_2['zenith_angle'])[0]

        t1 = list(zenith_1['bright_temp%s (K)' %f])[0]
        t2 = list(zenith_2['bright_temp%s (K)' %f])[0]
        t3 = list(zenith_3['bright_temp%s (K)' %f])[0]


        part_1 =  1 / ( sec(a2) - sec(a1) )

        a = t1 - t2
        b = t2 - t3

        part_2 = np.log( a / b )


        opacity = part_1 * part_2

        opac.append(opacity)


        # mean temp

        temp = (t1 * t3 - (t2**2) ) / (t1 + t3 - (2 * t2))

        #print(temp)

        mean_temp.append(temp)


    #zenith_1['opac_23'] = [opac[0]]
    #zenith_1['opac_31'] = [opac[1]]

    #print(opac, mean_temp)

    return opac, mean_temp


def tau_to_pwv(df, chan):


    if chan == 23:

        A = -59.697495
        B = 147.0750
        C = -1.8127

    elif chan == 31:

        A = -4369.135961
        B = 936.3715
        C = -15.7715

    else:

        pass


    tau_ = df['tau%s' %chan]

    pwv = (A * (tau_ ** 2)) + (B * tau_) + C

    df['PWV (mm)_%s' %chan] = pwv


    return df





def LineFit(x_data, y_data):

    slope, intercept, r, p, std_err = stats.linregress(x_data, y_data)

    #plt.plot(x_data, y_data)
    #plt.show()

    return slope, intercept




#Retrieving the antenna temperature using eq 5 from the technical documentation
def RetrieveAntTemp(data_frame, Th, chan):

    #return ((Pa - Pc)*(Ph - Pc)/(Th - Tc)) + Tc


    Pa = data_frame['%sG Sky Voltage (V)' %chan]
    Ph = data_frame['%sG Noise Source Voltage (V)' %chan]
    Pc = data_frame["%sG Load Voltage (V)" %chan]
    Tc = data_frame['%sG Plate Temperature (K)' %chan]

    a_g = (Pa - Pc)/(Ph - Pc)
    sec =  Th - Tc
    ant_temp = a_g * sec + Tc

    data_frame['ant_temp%s (K)' %chan] = ant_temp


    return data_frame




#Retrieving the antenna temperature using eq 5 from the technical documentation
def RetrieveAntTemp_final(data_frame, Th_chan):

    #return ((Pa - Pc)*(Ph - Pc)/(Th - Tc)) + Tc

    channels = [23, 31]
    Th_ = Th_chan


    chan_ant = []

    for chan, Th in zip(channels, Th_):

        Pa = data_frame['%sG Sky Voltage (V)' %chan]
        Ph = data_frame['%sG Noise Source Voltage (V)' %chan]
        Pc = data_frame["%sG Load Voltage (V)" %chan]
        Tc = data_frame['%sG Plate Temperature (K)' %chan]



        a_g = (Pa - Pc)/(Ph - Pc)
        sec =  Th - Tc
        ant_temp = a_g * sec + Tc

        chan_ant.append(ant_temp)

    data_frame['bright_temp23 (K)'] = chan_ant[0]
    data_frame['bright_temp31 (K)'] = chan_ant[1]


    zen_col = data_frame.pop('zenith_angle')
    data_frame.insert(18, 'zenith_angle', zen_col)

    print(data_frame)

    #zenith_0 = data_frame[data_frame['zenith_angle'].between(-1, 3)]

    #zenith_1 = data_frame[data_frame['zenith_angle'].between(28, 32)]

    #zenith_2 = data_frame[data_frame['zenith_angle'].between(38, 42)]

    #T_0 = zenith_0['bright_temp23 (K)']

    #zenith_1 = data_frame.loc[data_frame['bright_temp23 (K)'] > float(T_0) ].iloc[:1]

    #T_1 = zenith_1['bright_temp23 (K)']

    #zenith_2 = data_frame.loc[data_frame['bright_temp23 (K)'] > float(T_1) ].iloc[:1]

    #print(data_frame)

    #opacity, mean_temp = opacity_mean_temp(zenith_0, zenith_1, zenith_2)
    opacity = 1

    return data_frame, opacity




#Calculate airmass
def Airmass(A):

    A['airmass'] = 1 / np.sin( np.pi * A['Elevation angle (Deg)'] / 180 )

    #print(A)

    return A



def TipCurveCalibrate(dataline, chan):


    Th = 700
    error = 1000
    accuracy = 0.001

    while abs(error) > accuracy:

        retrieved_dataframe = RetrieveAntTemp(dataline, Th, chan)

        airmass = retrieved_dataframe['airmass']
        temp = retrieved_dataframe['ant_temp%s (K)' %chan]

        slope, intercept = LineFit(airmass, temp)


        error = intercept - 2.75



        Th = Th + error/10


    return Th











def read_in_and_process(file):

    df = pd.read_csv(file ,header=None)

    try:

        column_list = ["Sample Number",
        "Date & Time",
        "Elevation angle (Deg)",
        "31G Load Voltage (V)",
        "31G Noise Source Voltage (V)",
        "31G Sky Voltage (V)",
        "23G Load Voltage (V)",
        "23G Noise Source Voltage (V)",
        "23G Sky Voltage (V)",
        "23G Noise Source Temperature (K)",
        "23G Plate Temperature (K)",
        "31G Noise Source Temperature (K)",
        "31G Plate Temperature (K)",
        "Outdoor Temperature (K)",
        "Outdoor Relative Humidity (%)",
        "Barometric Pressure (mbar)",
        "Radometer internal temperature (C)"]

        df.columns = column_list

        #print(df)

    except:


        column_list = ["Sample Number",
        "Date & Time",
        "Elevation angle (Deg)",
        "31G Load Voltage (V)",
        "31G Noise Source Voltage (V)",
        "31G Sky Voltage (V)",
        "23G Load Voltage (V)",
        "23G Noise Source Voltage (V)",
        "23G Sky Voltage (V)",
        "23G Noise Source Temperature (K)",
        "23G Plate Temperature (K)",
        "31G Noise Source Temperature (K)",
        "31G Plate Temperature (K)",
        "Outdoor Temperature (K)",
        "Outdoor Relative Humidity (%)",
        "Barometric Pressure (mbar)"]

        df.columns = column_list

        #print(df)

    df = df.set_index("Date & Time")
    df = df.sort_index()



    return df


############################################# here the processing starts ################


try:

    parser=argparse.ArgumentParser()
    parser.add_argument('-p','--path', help='path to directory with str', type=str)
    args = parser.parse_args()
    path_to_files = args.path

except:

    pass # if no command path is added, then pass and path is None, empty. The default path will be used.


try:

    rad22_filenames = glob(path_to_files) ## reading in all the file names

except:

    ## This is the defualt path to file when empty arg path is provided on command

    #path_to_files = 'T240312.csv'   #path to one file, file sent by Reuben, uncomment for one file and comment next line
    path_to_files = "/Users/thecuriousvambo/22_rad/all_unprocces_22/*.csv" # path to all files
    rad22_filenames = glob(path_to_files) ## reading in all the file names


all_data_frames = []   ## empty list to store dataframes

for names in rad22_filenames:    ## looping through all the files
    all_data_frames.append(read_in_and_process(names))   ### appending all the empty files

all_data = pd.concat(all_data_frames, axis=0)    ## combined all the datafiles
all_data = all_data[all_data.index.notnull()]  ## removing nan values

df_to_work = all_data.sort_index()   ## sorting the dataframe according to date

group_ids = (df_to_work["Elevation angle (Deg)"] > (df_to_work["Elevation angle (Deg)"].shift() )).cumsum()  ## grouping by elevation for a fulltip curve
grouped = df_to_work.groupby(group_ids)   ## group each single tipcurve


all_zenith = []    ##  empty list for storing zenith direction dataframes

for k,g in grouped:   ## looping through each full tipcurve dataframe

     if len(g) == 8:   ## only considering a dataframe with a full tipcurve of 8 measurments

         print("----- Full tip curve dataframe ------")

         full_tip = Airmass(g)   ## calculating the airmass by calling airmass function

         full_tip['zenith_angle'] = 90 - full_tip['Elevation angle (Deg)']  ## calculating of the zenith angle

         #full_tip = full_tip[(full_tip['Elevation angle (Deg)'].between(88, 91))| (full_tip['Elevation angle (Deg)'].between(27, 29)) | (full_tip['Elevation angle (Deg)'].between(29, 32))]

         tip_linear_23 = TipCurveCalibrate(full_tip, 23)   ## estimating the antenna temp for 23 channel

         tip_linear_31 = TipCurveCalibrate(full_tip, 31)    ## estimating the antenna temp for 31 channel

         brightness_temp, zenith_bright = RetrieveAntTemp_final(full_tip, [tip_linear_23, tip_linear_31])  ## calculating the brightness temp


        # zenith_vals = brightness_temp[brightness_temp['zenith_angle'].between(-1, 2)]


         #all_zenith.append(zenith_vals)


#all_zen = pd.concat(all_zenith, axis=0)
#all_zen = all_zen[all_zen.index.notnull()]
#all_zen.to_csv('all_zen_22.csv')
