# -*- coding: utf-8 -*-
"""
Created on Sun Jun 15 19:41:46 2025

@author: jrobs

Initial Data Wrangling for the Atmospheric Measurements of Greenhouse Gases in the High Artic

Flask (and Insitu once I have it) data from ALERT, BARROW and SUMMIT research stations

1. Read in text files and concatenate all together, add "type" and "datetime" column for sorting
** Need to make the DataFrames better, each type of compound gets it's own column and nans 
fill the empty spaces

2. Plot the data from the 3 reasearch stations

3. Add the NEEM scenarios onto the plots
"""

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib qt
import seaborn as sns
import os
import glob

#switch working directory to the USRA_2025
cwd = os.getcwd()
wd = 'C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Code'
if cwd != wd:
    os.chdir(wd)
#%% grab all of the data

co2_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\co2\\*.txt"
co2_files = glob.glob(co2_path)

ch4_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\ch4\\*.txt"
ch4_files = glob.glob(ch4_path)

sf6_path= "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\sf6\\*.txt"
sf6_files = glob.glob(sf6_path)

ccl4_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\ccl4\\*.txt"
ccl4_files= glob.glob(ccl4_path)

ch3ccl3_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\ch3ccl3\\*.txt"
ch3ccl3_files = glob.glob(ch3ccl3_path)

cfc11_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\cfc11\\*.txt"
cfc11_files = glob.glob(cfc11_path)

cfc12_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\cfc12\\*.txt"
cfc12_files = glob.glob(cfc12_path)

cfc113_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\cfc113\\*.txt"
cfc113_files = glob.glob(cfc113_path)

#%% Plot the Co2 data together

for file in co2_files:
    data = pd.read_csv(file)
    

#%% concatenate all the alert data files together
for i,d in enumerate(alt_data_files):
    print(data_types[i], d)
    data = np.loadtxt(d, usecols=(1,2,3))
    df = pd.DataFrame(data, columns = ["year", "month", "value"])
    df["type"] = data_types[i]
    df["day"] = 1 #need for datetime object
    #make datetime column for plotting
    df["datetime"] = pd.to_datetime(df[["year", "month","day"]])
    
    if i == 0:
        alt_data_raw = df
    else:
        alt_data_raw = pd.concat([alt_data_raw, df])

#%% concatenate all the barrow files together
for i,d in enumerate(brw_data_files):
    print(data_types[i], d)
    data = np.loadtxt(d, usecols=(1,2,3))
    df = pd.DataFrame(data, columns = ["year", "month", "value"])
    df["type"] = data_types[i]
    df["day"] = 1 #need for datetime object
    #make datetime column for plotting
    df["datetime"] = pd.to_datetime(df[["year", "month","day"]])
    
    if i == 0:
        brw_data_raw = df
    else:
        brw_data_raw = pd.concat([brw_data_raw, df])
        
#%% concatenate all the summit data files together
for i,d in enumerate(sum_data_files):
    print(data_types[i], d)
    data = np.loadtxt(d, usecols=(1,2,3))
    df = pd.DataFrame(data, columns = ["year", "month", "value"])
    df["type"] = data_types[i]
    df["day"] = 1 #need for datetime object
    #make datetime column for plotting
    df["datetime"] = pd.to_datetime(df[["year", "month","day"]])
    
    if i == 0:
        sum_data_raw = df
    else:
        sum_data_raw = pd.concat([sum_data_raw, df])
        
#%% Read in NEEM files
neem_filepath = 
#%% first pass at plotting
#plot the data
f, axs = plt.subplots(6, figsize = (3,12))
for i,s in enumerate(data_types):
    
    gg = alt_data_raw[data_raw["type"]== s]
    axs[i].scatter(x = gg["datetime"], y = gg["value"], s = 5)
    axs[i].set_title(s, fontsize = 16)
    plt.tight_layout()
    
#%% figure 3 (ch4, methane)
#plot the concentration of methane
f = plt.figure(3, clear = True, figsize = (12,3))
ch4_alt = alt_data_raw[alt_data_raw["type"] == "ch4"]
ch4_brw = brw_data_raw[brw_data_raw["type"]== "ch4"]
ch4_sum = sum_data_raw[sum_data_raw["type"]== "ch4"]

plt.plot(ch4_alt["datetime"], ch4_alt["value"], color = "r", label = "Alert")
plt.plot(ch4_brw["datetime"], ch4_brw["value"], color = "k", label = "Barrow")
plt.plot(ch4_sum["datetime"], ch4_sum["value"], color = "g", label = "Summit")

plt.legend()
plt.xlabel("Time")
plt.ylabel("$CH_4$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n Methane in the High Arctic")

#%% figure 4 (ch4c13, carbon fractionation)
f = plt.figure(4, clear = True, figsize = (12,3))
ch4c13_alt = alt_data_raw[alt_data_raw["type"] == "ch4c13"]
ch4c13_brw = brw_data_raw[brw_data_raw["type"]== "ch4c13"]
ch4c13_sum = sum_data_raw[sum_data_raw["type"]== "ch4c13"]

plt.plot(ch4c13_alt["datetime"], ch4c13_alt["value"], color = "r", label = "Alert")
plt.plot(ch4c13_brw["datetime"], ch4c13_brw["value"], color = "k", label = "Barrow")
plt.plot(ch4c13_sum["datetime"], ch4c13_sum["value"], color = "g", label = "Summit")

plt.legend()
plt.xlabel("Time")
plt.ylabel("$d^{13}CH_4$ (ppt)")
plt.title("Monthly Average of $^{13}C/^{12}C$ of \n Methane in the High Arctic")

#%% figure 5 (carbon monoxide, co)
f = plt.figure(5, clear = True, figsize = (12,3))
co_alt = alt_data_raw[alt_data_raw["type"] == "co"]
co_brw = brw_data_raw[brw_data_raw["type"]== "co"]
co_sum = sum_data_raw[sum_data_raw["type"]== "co"]

plt.plot(co_alt["datetime"], co_alt["value"], color = "r", label = "Alert")
plt.plot(co_brw["datetime"], co_brw["value"], color = "k", label = "Barrow")
plt.plot(co_sum["datetime"], co_sum["value"], color = "g", label = "Summit")

plt.legend()
plt.xlabel("Time")
plt.ylabel("$CO$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n Carbon Monoxide in the High Arctic")

#%% figure 6 (carbon dioxide, co2)
f = plt.figure(6, clear = True, figsize = (12,3))
co2_alt = alt_data_raw[alt_data_raw["type"] == "co2"]
co2_brw = brw_data_raw[brw_data_raw["type"]== "co2"]
co2_sum = sum_data_raw[sum_data_raw["type"]== "co2"]

plt.plot(co2_alt["datetime"], co2_alt["value"], color = "r", label = "Alert")
plt.plot(co2_brw["datetime"], co2_brw["value"], color = "k", label = "Barrow")
plt.plot(co2_sum["datetime"], co2_sum["value"], color = "g", label = "Summit")

plt.legend()
plt.xlabel("Time")
plt.ylabel("$CO_2$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n Carbon Dioxide in the High Arctic")

#%% figure 7 (molecular hydrogen, h2)
f = plt.figure(7, clear = True, figsize = (12,3))
h2_alt = alt_data_raw[alt_data_raw["type"] == "h2"]
h2_brw = brw_data_raw[brw_data_raw["type"]== "h2"]
h2_sum = sum_data_raw[sum_data_raw["type"]== "h2"]

plt.plot(h2_alt["datetime"], h2_alt["value"], color = "r", label = "Alert")
plt.plot(h2_brw["datetime"], h2_brw["value"], color = "k", label = "Barrow")
plt.plot(h2_sum["datetime"], h2_sum["value"], color = "g", label = "Summit")

plt.legend()
plt.xlabel("Time")
plt.ylabel("$H_2$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n Molecular Hydrogen in the High Arctic")

#%% figure 8 (dinitrogen monoxide, n2o)
f = plt.figure(8, clear = True, figsize = (12,3))
n2o_alt = alt_data_raw[alt_data_raw["type"] == "n2o"]
n2o_brw = brw_data_raw[brw_data_raw["type"]== "n2o"]
n2o_sum = sum_data_raw[sum_data_raw["type"]== "n2o"]

plt.plot(n2o_alt["datetime"], n2o_alt["value"], color = "r", label = "Alert")
plt.plot(n2o_brw["datetime"], n2o_brw["value"], color = "k", label = "Barrow")
plt.plot(n2o_sum["datetime"], n2o_sum["value"], color = "g", label = "Summit")

plt.legend()
plt.xlabel("Time")
plt.ylabel("$N_2 O$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n Dinitrogen Monoxide in the High Arctic")

#%% figure 9 (sulfur hexafluoride, sf6)
f = plt.figure(9, clear = True, figsize = (12,3))
sf6_alt = alt_data_raw[alt_data_raw["type"] == "sf6"]
sf6_brw = brw_data_raw[brw_data_raw["type"]== "sf6"]
sf6_sum = sum_data_raw[sum_data_raw["type"]== "sf6"]

plt.plot(sf6_alt["datetime"], sf6_alt["value"], color = "r", label = "Alert")
plt.plot(sf6_brw["datetime"], sf6_brw["value"], color = "k", label = "Barrow")
plt.plot(sf6_sum["datetime"], sf6_sum["value"], color = "g", label = "Summit")

plt.legend()
plt.xlabel("Time")
plt.ylabel("$SF_6$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n Sulfur Hexafluoride in the High Arctic")







        



    

        
    
            
            
            
            
        


