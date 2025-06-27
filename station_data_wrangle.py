# -*- coding: utf-8 -*-
"""
Created on Sun Jun 15 19:41:46 2025

@author: jrobs

    Update Atmospheric Scenarios for measured species at EastGRIP
--------------------------------------------------------------------------------

    Atospheric Species Measured: CO2, CH4, SF6, CCl4, CH3CCL3, CFC11, CFC12, CFC113

--- INPUT DATA:
    Recent data (post NEEM) are gathered from https://gml.noaa.gov/data/data.php in monthly resolution in text
    file format. Both insitu and flask data can ben used (only inistu data exists for halocompounds). 
    Summit (SUM) is thepreferred station for it's proximity to the EGRIP site, however both Alert (ALT) and 
    Barrow (BRW) can be used. The notebook contains code to calculate the offset/ratio to correct the 
    ALT and BRW data to Summit values.

    For EastGRIP we cut off at 2018.63 (july) and correct with Summit as a reference.

--- TIME COLUMN:
    Data from NEEM have dates in decimal format whereas NOAA gives year, month format, the code will add a 
    date column to the data using the decimal system from NEEM, for ease I have just put the correct month 
    ratio in a list and add the correct ratio to the year.

--- FUNCTIONS:
    There are 3 functions that can be applied to all datasets. They all work toward finding the correct ratio and 
    offset to correct data to Summit values. More information on each function is available in the docstrings below.

        make_data_file: organises the data along a time column, adds nans where data doesn't exist

        find_seasonal_cycle: subtracts the yearly average from each data point to find the mean seasonal cycle

        subtract_seasonal_cycle: subtracts the seasonal cycle from the monthly average data to remove the seasonal 
        cycle
        
--- CODE
    I have annotated the rest of the code where applicable and where anything out of the ordinary was done. The most
    annotations are in the CO2 section and subsequent sections are just similar versions of the CO2 workflow. ALl 
    sections use the same functions and produce the same plot. 
    
    A linear interpolation is applied in [HALOCOMPOUNDS] since there is missing data in the record.

--- CORRECTION:
    Once the data is in the correct format, we can take the linear regression of summit vs alert or barrow and use the 
    slope and intercept to correct the data. Some data will require only the offset, or only the ratio, so make sure
    to run the code that plots the corrected data before the averaging to make sure that it wasn't over corrected.

    *** For EastGRIP we correct to Summit, which is the closest, for Muller we would have to correct to Alert. 
    These sections would have to all be changed since I do the correction by hand instead of using a function

--- PLOTS:
    This code generates 6 plots for each species

        figure1: raw data at each station from 2000 to 2019 to find any problems with the data
    
        figure2: the mean seasonal cycle at each station with the actual values plotted as dots, check to see 
        if the is any stations that are glaringly different
    
        figure3: Monthly Average data with the seasonal cycle removed to make sure that the cycle was actually removed
    
        figure4: summit vs barrow and alert with the best fit lines
    
        figure5: Corrected data from 2009-2019 to check the Summit fit
    
    figure6: Averaged data from 1995-2015 with NEEM data to check the fit of the overlap

--- FINAL DATA:
        The final dataset (NEEM + new stuff) and all figures are saved in final_data_directory which would need to 
    be changed for any subsequent projects
'
"""

#import packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
%matplotlib qt 
#lets us have the interact-able figures
import seaborn as sns
import os
import glob
import scipy.stats as stats

#switch working directory to the USRA_2025, this is specifically because my computer hates me
cwd = os.getcwd()
wd = 'C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Code'
final_data_directory = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\Final\\co2\\"
if cwd != wd:
    os.chdir(wd)
    
date_addons = [0.04, 0.13, 0.21, 0.29, 0.38, 0.46, 0.54, 0.63, 0.71, 0.79, 0.88, 0.96]
#%% define functions
#make function to add date column
def add_dates(array, percents = date_addons):
    """"Function to add decimal dates to the NOAA data. Mean date is the middle of the month. Here the function pulls
    the correct ratio from a list by indexing for simplicity but the calculation is made as:
        
        decimal date = year + dayofyear/365.25 
        
        which is date_addons = [0.04, 0.13, 0.21, 0.29, 0.38, 0.46, 0.54, 0.63, 0.71, 0.79, 0.88, 0.96]"""
    array["date"] = 0
    for i, x in enumerate(percents):
        years_at_month =  array[array["month"]==i+1]["year"]
        percent = np.ones(len(years_at_month))*x
        array.loc[array["month"]==i+1, "date"] = (years_at_month + percent).values
    return array

def make_data_file(alt_data, brw_data, sum_data, time, month, year):
    """"Take data from 3 stations and arrange them on one time axis. The data must have columns "value" and 
    "date". The function will produce a 6 column array with date, year, month, alert, barrow and summit data. 
    When no data exists for a certain month at a station, the function will fill the spot in the dataframe with 
    a nan. 
    
    IMPORTANT** If using both insitu and flask data the data must be mixed prior to using this function, the function
    will only take 3 datasets as inputs.
    
    Inputs
    ----------------
    alt_data: pandas dataframe
        monthly average data from the alert station
    brw_data: pandas dataframe
        monthly average data from the barrow station
    sum_data: pandas dataframe
        monthly average data from the summit station
    time: numpy array
        decimal dates for the time period wanted, this will form the time axis for future plots.
    month: numpy array
        an array of months corresponding to the correct measurement. There should be one month for each row of the 
        dataset (ex. 2008.96 should have a month value of 12).
    year: numpy array
        array of years corresponding to the correct measurement. There should be one year for each row of the 
        dataset (ex. 2008.96 should have a year value of 2008)
        
    Returns
    ------------------
    a N X 6 dataframe where N is the length of the time period. 3 time columns and 3 station value columns will be
    provided."""
    
    #make data frame
    df = pd.DataFrame(time)
    df["year"] = year
    df["month"] = month
    df["alt"] = 0
    df["brw"] = 0
    df["sum"] = 0
    
    for i, t in enumerate(df["date"]):
        #grad the data at that specific date
        alt_at_t = alt_data[alt_data["date"] == t]["value"].values
        brw_at_t = brw_data[brw_data["date"] == t]["value"].values
        sum_at_t = sum_data[sum_data["date"] == t]["value"].values
        
        #check that the value at that specific date exists, else put a nan in the big data frame
        
        if alt_at_t.size > 0: #this will be true if there is a datapoint for this date
            df.loc[i,"alt"] = alt_at_t
        else:
            df.loc[i,"alt"] = np.nan #if the test fails put a nan instead
        
        if brw_at_t.size > 0:
            df.loc[i, "brw"] = brw_at_t
        else: 
            df.loc[i, "brw"] = np.nan

        if sum_at_t.size > 0:
            df.loc[i, "sum"] = sum_at_t
        else:
            df.loc[i, "sum"] = np.nan
    
    return df
def find_seasonal_cycle(df, stations = ["alt", "brw", "sum"]):
    """Finds the seasonal cycle for 3 stations in a dataset. 
    If using the make_data_file function, this function will susbtract the yearly average from the data and add 3 
    columns with the subsequent seasonal cycles. 
    
   If there is less than 12 months of data for a year (i.e there is data missing) the function will put nans in
   for that date in the dataframe. For our purposes we want complete seasonal cycles and missing data could 
   skew the data up or down when subtracting the yearly average. 
   
   Inputs
   ---------------
   df: pandas dataframe
       dataframe produced using the make_data_file function. A N x 6 dataframe with a date, year and month column
       along with monthly average data for alert, barrow and summit stations in the last 3 columns.
       
    stations: list, default: ['alt', 'brw', 'sum']
        the stations that the data was collected from. Also the names of the data columns in df
        
    Returns
    --------------
    The inputted dataframe with 3 extra columns containing the yearly anomalies of the data.
    """
    
    for s in stations:
        station_data = df[[s, "year", "month"]]
        df[s + "_anom"] = 0 #make anomaly column
        
        df_yearly = station_data.groupby('year').mean() #find the yearly mean for that year and station
    
        for y in df_yearly.index: #year
            value_at_year = station_data[station_data["year"]==y][s]
            if (value_at_year/value_at_year).sum() == 12: 
                #check to see if we have a full year. Nans are not counted as numbers so the test will fail
                df.loc[df["year"]==y, s+"_anom"] = value_at_year-df_yearly[s][y] 
                #subtract to make anomaly
            else: df.loc[df["year"]==y, s+"_anom"] = np.nan #make the year nans if the dataset is not complete
        
    return df

def clean_data(df, start_time = 2000, end_time = 2019):
    """ workflow to ease cleaning the data up for future functions. Calls add_dates to add a date column as well
   as croppping the data since we don't need the full record. Here we use 9 years of overlap since NEEM ends at 2008.96 """
    df = add_dates(df) #get dates in clean format
    #df = find_seasonal_cycle(df, plot = False) #get anomalies
    df = df[df["date"] > start_time] #grab all of the data that NEEM doesn't have with 10 years overlap
    df = df[df["date"] < end_time] #don't need all the way to 2023 for EGRIP project
    df = df.set_index(np.arange(0,len(df["value"]+1), 1))
    
    return df
    
def subtract_seasonal_cycle(df, df_seasonal_means, stations = ["alt", "brw", "sum"]):
    """
    Subtract the mean seasonal cycle from monthly averaged data. 

    Parameters
    ----------
    df : pandas dataframe
        Dataframe generated using make_data_file and find_seasonal_cycle functions.
        Has 3 columns with decimal dates, year and months as well as 3 columns of raw data from NOAA, and 3 columns
        with the yearly average subtracted to make the seasonal cycle.
    df_seasonal_means : pandas dataframe
        dataframe of the average seasonal cycle for each station. 
        take the monthly seasonal cycle, groupby month and take the mean for each station.
    stations : list, optional
        The stations that will be looped over. The default is ["alt", "brw", "sum"].

    Returns
    -------
    df : pandas dataframe
        Inputted dataframe with 3 added columns that contain the station data minus the average seasonal cycle at
        that station.

    """
    for s in stations:
        station_data = df[["year", "month", s]]
        df[s+"_yearly"] = 0
        for m in df_seasonal_means.index: #months
            cycle_at_month = df_seasonal_means[s+"_anom"][m]
            data_at_month = station_data[station_data["month"] == m]
            df.loc[df["month"]==m, s+"_yearly"] = data_at_month[s] - cycle_at_month
    return df
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


#%% useful stuff to keep things consistent
end_time = 2008.96 #last time point for NEEM data
stations = ["alt", "brw", "sum"] #weather station ids for reference

#which columns to grab for each type of data
#want the year, month and value (+std if available)
gg_columns_flask = (1,2,3) #the first column is the station_id and throws an error when read-in
gg_columns_insitu = (1,2,10,11) #this data is in a different format, grab the same columns + the standard deviations
hal_columns = (0,1,2,3,4) # for the halocompounds since those are also in a different form

#Make colormap for plots
cmap = mpl.colormaps['managua']
n_lines = 6

colors = cmap(np.linspace(0,1,n_lines))
#%% Get CO2 Data
co2_alt_flask = pd.DataFrame( np.loadtxt(co2_files[0], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
co2_alt_flask = clean_data(co2_alt_flask)


co2_brw_flask = pd.DataFrame( np.loadtxt(co2_files[1], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
co2_brw_flask = clean_data(co2_brw_flask)


co2_brw_insitu = pd.DataFrame( np.loadtxt(co2_files[2], usecols = gg_columns_insitu), columns = 
                             ["year", "month", "value", "std"])
co2_brw_insitu = clean_data(co2_brw_insitu)


co2_sum_flask = pd.DataFrame( np.loadtxt(co2_files[3], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
co2_sum_flask = clean_data(co2_sum_flask)


co2_NEEM = pd.DataFrame( np.loadtxt(co2_files[4], ), columns = 
                      ["date","value", "std"])

#average both barrow datasets together

co2_brw = pd.DataFrame((co2_brw_flask["value"] + co2_brw_insitu["value"])/2)
co2_brw["date"] = co2_brw_flask["date"].values
co2_brw.loc[34, "value"] = co2_brw_flask["value"][34]

#%% make a dataframe for all of the co2 data, plot the raw data
co2_data = make_data_file(co2_alt_flask, co2_brw, co2_sum_flask, time = co2_brw_flask["date"], year = co2_brw_flask["year"], 
month = co2_brw_flask["month"])       

plt.figure(1, clear = True)
plt.plot(co2_data["date"], co2_data["alt"], label = "ALERT", color = colors[0])
plt.plot(co2_data["date"], co2_data["brw"], label = "BARROW", color = colors[2])
plt.plot(co2_data["date"], co2_data["sum"], label = "SUMMIT", color = colors[4])

plt.xlabel("time")
plt.ylabel("$CO_2$ concentration (ppm)")
plt.title("Atmospheric Concentration of $CO_2$ in the High Arctic")
plt.legend(title = "--Station--")

plt.savefig(final_data_directory + "figure1.png")
#%% seasonal cycle of co2
co2_seasonal = find_seasonal_cycle(co2_data)
co2_seasonal_means = co2_seasonal.groupby("month").mean()
plt.figure(2, clear = True)
#plot the mean seasonal cycle
plt.plot(co2_seasonal_means.index, co2_seasonal_means["alt_anom"], color = colors[0], label = "alert")
plt.plot(co2_seasonal_means.index, co2_seasonal_means["brw_anom"], color = colors[2], label = "barrow")
plt.plot(co2_seasonal_means.index, co2_seasonal_means["sum_anom"], color = colors[4], label = "summit")
#plot each individual point
plt.scatter(co2_seasonal["month"], co2_seasonal["alt_anom"], color = colors[0], alpha = 0.5)
plt.scatter(co2_seasonal["month"], co2_seasonal["brw_anom"], color = colors[2], alpha = 0.5)
plt.scatter(co2_seasonal["month"], co2_seasonal["sum_anom"], color = colors[4], alpha = 0.5)

plt.xlabel("month")
plt.ylabel("yearly anomaly (ppm)")
plt.title("Seasonal Cycle of $CO_2$ at 3 stations in the High Arctic")
plt.legend(loc = "lower left", title = "--Station--")

plt.savefig(final_data_directory + "figure2.png")
#%% subtract the seasonal cycle from each station
co2_data = subtract_seasonal_cycle(co2_data, co2_seasonal_means)
plt.figure(3, clear = True)
plt.plot(co2_data["date"], co2_data["alt_yearly"], color = colors[0], label = "alert")
plt.plot(co2_data["date"], co2_data["brw_yearly"], color = colors[2], label = "brw")
plt.plot(co2_data["date"], co2_data["sum_yearly"], color = colors[4], label = "summit")

plt.xlabel("time")
plt.ylabel("concentration (ppm)")
plt.title("Concentration of Atmospheric $CO_2$ in the High Arctic \n with the Seasonal Cycle Removed")
plt.legend(title = "--Station--")

plt.savefig(final_data_directory + "figure3.png")
#%% plot the ratio of each station against summit to find the ratio and offset
co2_data_nona = co2_data.dropna()
alt_regress = stats.linregress(co2_data_nona["sum_yearly"], co2_data_nona["alt_yearly"])
brw_regress = stats.linregress(co2_data_nona["sum_yearly"], co2_data_nona["brw_yearly"])

alt_ratio = alt_regress.slope
brw_ratio = brw_regress.slope
alt_offset = alt_regress.intercept #I'm skeptical of this offset so I didn't include it when I corrected the data
brw_offset = brw_regress.intercept

plt.figure(4, clear = True)
plt.scatter(co2_data["sum_yearly"], co2_data["alt_yearly"], label = "alert", color = colors[0], alpha = 0.7)
plt.scatter(co2_data["sum_yearly"], co2_data["brw_yearly"], label = "barrow", color = colors[2], alpha = 0.7)

alert_line = co2_data["sum_yearly"]*alt_regress.slope + alt_regress.intercept
brw_line = co2_data["sum_yearly"]*brw_regress.slope + brw_regress.intercept

plt.plot(co2_data["sum_yearly"], alert_line, color = colors[0])
plt.plot(co2_data["sum_yearly"], brw_line, color = colors[2])

plt.xlabel("summit $CO_2$ (ppm)")
plt.ylabel("other station's concentration (ppm)")
plt.legend()

plt.annotate(f'ALERT : y = {round(alt_ratio,4)}x + {round(alt_offset,4)}', xy = (400,380), color = "k")
plt.annotate(f'BARROW : y = {round(brw_ratio,4)}x + {round(brw_offset,4)}', xy = (400,378), color = "k")

plt.title("Regression of Summit Against Barrow and Alert")
plt.savefig(final_data_directory+"figure4.png")

#%% Apply Offset and Ratio and Merge with NEEM

co2_egrip = co2_data[["date", "alt", "brw", "sum"]].copy()
co2_egrip["alt"] = co2_data["alt"]*alt_ratio #- alt_offset
co2_egrip["brw"] = co2_data["brw"]*brw_ratio - brw_offset

co2_egrip["value"] = 0
co2_egrip["std"] = 0 
for i in co2_egrip.index:
    data_at_i = co2_egrip.loc[i][stations]
    
    co2_egrip.loc[i,"value"] = data_at_i.mean()
    co2_egrip.loc[i, "std"] = data_at_i.std()

co2_overlap = co2_egrip[co2_egrip["date"] <= 2008.96] #make another dataframe to look at the overlapping 9 years
co2_egrip = co2_egrip[co2_egrip["date"] > 2008.96] #keep this for the final dataset
co2_egrip_final = pd.concat([co2_NEEM, co2_egrip[["date", "value", "std"]]]) #make the final dataset
co2_egrip_final.to_csv(final_data_directory+"SCENARIO_EGRIP18_CO2.csv")

#plot the new data
plt.figure(5, clear = True)

plt.plot(co2_egrip["date"], co2_egrip["alt"], color = colors[0], label = "alert")
plt.plot(co2_egrip["date"], co2_egrip["brw"], color = colors[2], label = "barrow")
plt.plot(co2_egrip["date"], co2_egrip["sum"], color = colors[4], label = "summit")
plt.plot(co2_NEEM[co2_NEEM["date"] > 2000]["date"], co2_NEEM[co2_NEEM["date"] > 2000]["value"], 
         color = colors[5], label = "NEEM")

plt.legend(title = "--station--")
plt.savefig(final_data_directory+"figure5.png")
#%% plot the final dataset
plt.figure(6, clear = True)
plt.plot(co2_egrip_final["date"], co2_egrip_final["value"], color = colors[3])
plt.fill_between(co2_egrip_final["date"], co2_egrip_final["value"] - co2_egrip_final["std"], 
                 co2_egrip_final["value"] + co2_egrip_final["std"], color = colors[3], label = "EGRIP", alpha = 0.5)
plt.xlim(1995, 2015)
plt.plot(co2_NEEM["date"], co2_NEEM["value"], color = colors[5])
plt.fill_between(co2_NEEM["date"], co2_NEEM["value"] - co2_NEEM["std"], 
                 co2_NEEM["value"] + co2_NEEM["std"], color = colors[5], label = "NEEM", alpha = 0.5)

plt.plot(co2_overlap["date"], co2_overlap["value"], color = colors[1])
plt.fill_between(co2_overlap["date"], co2_overlap["value"] - co2_overlap["std"], 
                 co2_overlap["value"] + co2_overlap["std"], color = colors[1], label = "EGRIP", alpha = 0.5)
plt.savefig(final_data_directory+"figure6.png")
#%% ########## CH4 DATA ##################
ch4_alt_flask = pd.DataFrame( np.loadtxt(ch4_files[0], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
ch4_alt_flask = clean_data(ch4_alt_flask)


ch4_brw_flask = pd.DataFrame( np.loadtxt(ch4_files[1], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
ch4_brw_flask = clean_data(ch4_brw_flask)


ch4_brw_insitu = pd.DataFrame( np.loadtxt(ch4_files[2], usecols = gg_columns_insitu), columns = 
                             ["year", "month", "value", "std"])
ch4_brw_insitu = clean_data(ch4_brw_insitu)


ch4_sum_flask = pd.DataFrame( np.loadtxt(ch4_files[3], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])

ch4_sum_flask = clean_data(ch4_sum_flask)


ch4_NEEM = pd.DataFrame( np.loadtxt(ch4_files[4], ), columns = 
                             ["date","value", "std"])

#plot the co2 data
plt.figure(2, clear = True)

plt.plot(ch4_alt_flask["date"],ch4_alt_flask["value"], 
         label = "ALT (flask)", color = colors[0])
plt.plot(ch4_brw_flask["date"],ch4_brw_flask["value"], 
         label = "BRW (flask)", color = colors[1])
plt.plot(ch4_sum_flask["date"], ch4_sum_flask["value"], 
         label = "SUM (flask)", color = colors[2])
plt.plot(ch4_brw_insitu["date"], ch4_brw_insitu["value"], 
         label = "BRW (insitu)", color = colors[3])
plt.plot(ch4_NEEM["date"], ch4_NEEM["value"],
         label = "NEEM", color = colors[4])
plt.legend()

plt.xlabel("Time")
plt.ylabel("$CH_4$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n Methane in the High Arctic")


#%% plot ch4 seasonal cycle
plt.figure(10, clear = True)
vmin = 2008.96
vmax = 2025
h = plt.scatter(ch4_alt_flask["month"], ch4_alt_flask["anom"], c = ch4_alt_flask["year"], cmap = cmap, marker = "o" ,
            label = "alert", vmin = vmin, vmax = vmax)
plt.scatter(ch4_brw_flask["month"], ch4_brw_flask["anom"], c = ch4_brw_flask["year"], cmap = cmap, marker = "v" ,
            label = "barrow", vmin = vmin, vmax = vmax)
plt.scatter(ch4_sum_flask["month"], ch4_sum_flask["anom"], c = ch4_sum_flask["year"], cmap = cmap, marker = "*" ,
            label = "summit", vmin = vmin, vmax = vmax)
plt.scatter(ch4_brw_insitu["month"], ch4_brw_insitu["anom"], c = ch4_brw_insitu["year"], cmap = cmap, marker = "s" ,
            label = "barrow (insitu", vmin = vmin, vmax = vmax)
cbar = plt.colorbar(h)
plt.xlabel("month")
plt.ylabel("anomaly (ppm)")
plt.title("Seasonal Cycle of ch4")
cbar.set_label("year")
plt.legend()
#%% Get Sf6 Data
sf6_alt_flask = pd.DataFrame( np.loadtxt(sf6_files[2], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
sf6_alt_flask = clean_data(sf6_alt_flask)


sf6_brw_flask = pd.DataFrame( np.loadtxt(sf6_files[3], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
sf6_brw_flask = clean_data(sf6_brw_flask)


sf6_brw_insitu = pd.DataFrame( np.loadtxt(sf6_files[0], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
sf6_brw_insitu = clean_data(sf6_brw_insitu)


sf6_sum_flask = pd.DataFrame( np.loadtxt(sf6_files[4], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
sf6_sum_flask = clean_data(sf6_sum_flask)


sf6_sum_insitu = pd.DataFrame( np.loadtxt(sf6_files[5], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
sf6_sum_insitu = clean_data(sf6_sum_insitu)


sf6_NEEM = pd.DataFrame( np.loadtxt(sf6_files[1]), columns = 
                             ["date","value", "std"])

#plot the co2 data
plt.figure(3, clear = True)

plt.plot(sf6_alt_flask["date"],sf6_alt_flask["value"], 
         label = "ALT (flask)", color = colors[0])
plt.plot(sf6_brw_flask["date"],sf6_brw_flask["value"], 
         label = "BRW (flask)", color = colors[1])
plt.plot(sf6_sum_flask["date"], sf6_sum_flask["value"], 
         label = "SUM (flask)", color = colors[2])
plt.plot(sf6_brw_insitu["date"], sf6_brw_insitu["value"], 
         label = "BRW (insitu)", color = colors[3])
plt.plot(sf6_NEEM["date"], sf6_NEEM["value"],
         label = "NEEM", color = colors[4])
plt.plot(sf6_sum_insitu["date"], sf6_sum_insitu["value"],
         label = "SUMMIT (insitu)", color = colors[5])
plt.legend()

plt.xlabel("Time")
plt.ylabel("$SF_6$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n $SF_6$ in the High Arctic")
#%% seasonal cycle of sf6
plt.figure(11, clear = True)
vmin = 2008.96
vmax = 2025
h = plt.scatter(sf6_alt_flask["month"], sf6_alt_flask["anom"], c = sf6_alt_flask["year"], cmap = cmap, marker = "o" ,
            label = "alert", vmin = vmin, vmax = vmax)
plt.scatter(sf6_brw_flask["month"], sf6_brw_flask["anom"], c = sf6_brw_flask["year"], cmap = cmap, marker = "v" ,
            label = "barrow", vmin = vmin, vmax = vmax)
plt.scatter(sf6_sum_flask["month"], sf6_sum_flask["anom"], c = sf6_sum_flask["year"], cmap = cmap, marker = "*" ,
            label = "summit", vmin = vmin, vmax = vmax)
plt.scatter(sf6_brw_insitu["month"], sf6_brw_insitu["anom"], c = sf6_brw_insitu["year"], cmap = cmap, marker = "s" ,
            label = "barrow (insitu", vmin = vmin, vmax = vmax)
plt.scatter(sf6_sum_insitu["month"], sf6_sum_insitu["anom"], c = sf6_sum_insitu["year"], cmap = cmap, marker = "^" ,
            label = "summit (insitu", vmin = vmin, vmax = vmax)
cbar = plt.colorbar(h)
plt.xlabel("month")
plt.ylabel("anomaly (ppm)")
plt.title("Seasonal Cycle of sf6")
cbar.set_label("year")
plt.legend()
#%% CCL4

ccl4_brw_insitu = pd.DataFrame( np.loadtxt(ccl4_files[0], usecols= hal_columns), 
                               columns = ["year", "month", "value", "unc", "std"])
ccl4_brw_insitu = clean_data(ccl4_brw_insitu)

ccl4_sum_insitu = pd.DataFrame( np.loadtxt(ccl4_files[2], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
ccl4_sum_insitu = clean_data(ccl4_sum_insitu)


ccl4_NEEM = pd.DataFrame( np.loadtxt(ccl4_files[1], ), columns = 
                             ["date","value", "std"])

#plot the co2 data
plt.figure(4, clear = True)

plt.plot(ccl4_brw_insitu["date"], ccl4_brw_insitu["value"], 
         label = "BRW (insitu)", color = colors[3])
plt.plot(ccl4_NEEM["date"], ccl4_NEEM["value"],
         label = "NEEM", color = colors[4])
plt.plot(ccl4_sum_insitu["date"], ccl4_sum_insitu["value"],
         label = "SUMMIT (insitu)", color = colors[5])
plt.legend()

plt.xlabel("Time")
plt.ylabel("$CCL_4$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n $CCL_4$ in the High Arctic")

#%% seasonal cycle of ccl4
plt.figure(12, clear = True)
vmin = 2008.96
vmax = 2025
h = plt.scatter(ccl4_brw_insitu["month"], ccl4_brw_insitu["anom"], c = ccl4_brw_insitu["year"], cmap = cmap, marker = "s" ,
            label = "barrow (insitu", vmin = vmin, vmax = vmax)
plt.scatter(ccl4_sum_insitu["month"], ccl4_sum_insitu["anom"], c = ccl4_sum_insitu["year"], cmap = cmap, marker = "^" ,
            label = "summit (insitu", vmin = vmin, vmax = vmax)
cbar = plt.colorbar(h)
plt.xlabel("month")
plt.ylabel("anomaly (ppm)")
plt.title("Seasonal Cycle of ccl4")
cbar.set_label("year")
plt.legend()
#%% CH3CCL3

ch3ccl3_brw_insitu = pd.DataFrame( np.loadtxt(ch3ccl3_files[0], usecols= hal_columns), 
                               columns = ["year", "month", "value", "unc", "std"])
ch3ccl3_brw_insitu = clean_data(ch3ccl3_brw_insitu)

ch3ccl3_brw_insitu1 = pd.DataFrame( np.loadtxt(ch3ccl3_files[1], usecols= hal_columns), 
                               columns = ["year", "month", "value", "unc", "std"])
ch3ccl3_brw_insitu1 = clean_data(ch3ccl3_brw_insitu1)

ch3ccl3_brw_insitu = pd.concat([ch3ccl3_brw_insitu, ch3ccl3_brw_insitu1], 
                               ignore_index= True)

ch3ccl3_sum_insitu = pd.DataFrame( np.loadtxt(ch3ccl3_files[2], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
ch3ccl3_sum_insitu = clean_data(ch3ccl3_sum_insitu)


ch3ccl3_NEEM = pd.DataFrame( np.loadtxt(ch3ccl3_files[3], ), columns = 
                             ["date","value", "std"])

#plot the co2 data
plt.figure(5, clear = True)

plt.plot(ch3ccl3_brw_insitu["date"], ch3ccl3_brw_insitu["value"], 
         label = "BRW (insitu)", color = colors[3])
plt.plot(ch3ccl3_NEEM["date"], ch3ccl3_NEEM["value"],
         label = "NEEM", color = colors[4])
plt.plot(ch3ccl3_sum_insitu["date"], ch3ccl3_sum_insitu["value"],
         label = "SUMMIT (insitu)", color = colors[5])
plt.legend()

plt.xlabel("Time")
plt.ylabel("$CH_3 CCl_3$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n $CH_3 CCl_3$ in the High Arctic")

#%% seasonal cycle of ch3ccl3

plt.figure(13, clear = True)
vmin = 2008.96
vmax = 2025
h = plt.scatter(ch3ccl3_brw_insitu["month"], ch3ccl3_brw_insitu["anom"], c = ch3ccl3_brw_insitu["year"], cmap = cmap, marker = "s" ,
            label = "barrow (insitu", vmin = vmin, vmax = vmax)
plt.scatter(ch3ccl3_sum_insitu["month"], ch3ccl3_sum_insitu["anom"], c = ch3ccl3_sum_insitu["year"], cmap = cmap, marker = "^" ,
            label = "summit (insitu", vmin = vmin, vmax = vmax)
cbar = plt.colorbar(h)
plt.xlabel("month")
plt.ylabel("anomaly (ppm)")
plt.title("Seasonal Cycle of ch3ccl3")
cbar.set_label("year")
plt.legend()

#%% CFC11
cfc11_brw_insitu = pd.DataFrame(np.loadtxt(cfc11_files[0], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc11_brw_insitu1 = pd.DataFrame(np.loadtxt(cfc11_files[1], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc11_brw_insitu = pd.concat([cfc11_brw_insitu, cfc11_brw_insitu1], ignore_index = True)

cfc11_sum_insitu = pd.DataFrame(np.loadtxt(cfc11_files[2], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc11_sum_insitu = clean_data(cfc11_sum_insitu)
cfc11_brw_insitu = clean_data(cfc11_brw_insitu)

cfc11_NEEM = pd.DataFrame(np.loadtxt(cfc11_files[4], usecols = (0,1,2)), columns = 
                                ["date", "value", "std"])

#plot the co2 data
plt.figure(6, clear = True)

plt.plot(cfc11_brw_insitu["date"], cfc11_brw_insitu["value"], 
         label = "BRW (insitu)", color = colors[3])
plt.plot(cfc11_NEEM["date"], cfc11_NEEM["value"],
         label = "NEEM", color = colors[4])
plt.plot(cfc11_sum_insitu["date"], cfc11_sum_insitu["value"],
         label = "SUMMIT (insitu)", color = colors[5])
plt.legend()

plt.xlabel("Time")
plt.ylabel("$CFC11$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n $CFC11$ in the High Arctic")

#%% seasonal cycle of cfc11
plt.figure(14, clear = True)
vmin = 2008.96
vmax = 2025
h = plt.scatter(cfc11_brw_insitu["month"], cfc11_brw_insitu["anom"], c = cfc11_brw_insitu["year"], cmap = cmap, marker = "s" ,
            label = "barrow (insitu", vmin = vmin, vmax = vmax)
plt.scatter(cfc11_sum_insitu["month"], cfc11_sum_insitu["anom"], c = cfc11_sum_insitu["year"], cmap = cmap, marker = "^" ,
            label = "summit (insitu", vmin = vmin, vmax = vmax)
cbar = plt.colorbar(h)
plt.xlabel("month")
plt.ylabel("anomaly (ppm)")
plt.title("Seasonal Cycle of cfc11")
cbar.set_label("year")
plt.legend()
#%% CFC12
cfc12_brw_insitu = pd.DataFrame(np.loadtxt(cfc12_files[0], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc12_brw_insitu1 = pd.DataFrame(np.loadtxt(cfc12_files[1], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc12_brw_insitu = pd.concat([cfc12_brw_insitu, cfc12_brw_insitu1], ignore_index = True)

cfc12_sum_insitu = pd.DataFrame(np.loadtxt(cfc12_files[2], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc12_sum_insitu = clean_data(cfc12_sum_insitu)
cfc12_brw_insitu = clean_data(cfc12_brw_insitu)

cfc12_NEEM = pd.DataFrame(np.loadtxt(cfc12_files[3], usecols = (0,1,2)), columns = 
                                ["date", "value", "std"])

#plot the co2 data
plt.figure(7, clear = True)

plt.plot(cfc12_brw_insitu["date"], cfc12_brw_insitu["value"], 
         label = "BRW (insitu)", color = colors[3])
plt.plot(cfc12_NEEM["date"], cfc12_NEEM["value"],
         label = "NEEM", color = colors[4])
plt.plot(cfc12_sum_insitu["date"], cfc12_sum_insitu["value"],
         label = "SUMMIT (insitu)", color = colors[5])
plt.legend()

plt.xlabel("Time")
plt.ylabel("$CFC12$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n $CFC12$ in the High Arctic")

#%% seasonal cycle of cfc12
plt.figure(15, clear = True)
vmin = 2008.96
vmax = 2025
h = plt.scatter(cfc12_brw_insitu["month"], cfc12_brw_insitu["anom"], c = cfc12_brw_insitu["year"], cmap = cmap, marker = "s" ,
            label = "barrow (insitu", vmin = vmin, vmax = vmax)
plt.scatter(cfc12_sum_insitu["month"], cfc12_sum_insitu["anom"], c = cfc12_sum_insitu["year"], cmap = cmap, marker = "^" ,
            label = "summit (insitu", vmin = vmin, vmax = vmax)
cbar = plt.colorbar(h)
plt.xlabel("month")
plt.ylabel("anomaly (ppm)")
plt.title("Seasonal Cycle of cfc12")
cbar.set_label("year")
plt.legend()
#%% CFC113
cfc113_brw_insitu = pd.DataFrame(np.loadtxt(cfc113_files[0], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])

cfc113_sum_insitu = pd.DataFrame(np.loadtxt(cfc113_files[1], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc113_sum_insitu = clean_data(cfc113_sum_insitu)
cfc113_brw_insitu = clean_data(cfc113_brw_insitu)

cfc113_NEEM = pd.DataFrame(np.loadtxt(cfc113_files[2], usecols = (0,1,2)), columns = 
                                ["date", "value", "std"])

#plot the co2 data
plt.figure(8, clear = True)

plt.plot(cfc113_brw_insitu["date"], cfc113_brw_insitu["value"], 
         label = "BRW (insitu)", color = colors[3])
plt.plot(cfc113_NEEM["date"], cfc113_NEEM["value"],
         label = "NEEM", color = colors[4])
plt.plot(cfc113_sum_insitu["date"], cfc113_sum_insitu["value"],
         label = "SUMMIT (insitu)", color = colors[5])
plt.legend()

plt.xlabel("Time")
plt.ylabel("$CFC113$ concentration (ppt)")
plt.title("Monthly Average Concentration of \n $CFC113$ in the High Arctic")

#%% seasonal cycle of cfc113

plt.figure(16, clear = True)
vmin = 2008.96
vmax = 2025
h = plt.scatter(cfc113_brw_insitu["month"], cfc113_brw_insitu["anom"], c = cfc113_brw_insitu["year"], cmap = cmap, marker = "s" ,
            label = "barrow (insitu", vmin = vmin, vmax = vmax)
plt.scatter(cfc113_sum_insitu["month"], cfc113_sum_insitu["anom"], c = cfc113_sum_insitu["year"], cmap = cmap, marker = "^" ,
            label = "summit (insitu", vmin = vmin, vmax = vmax)
cbar = plt.colorbar(h)
plt.xlabel("month")
plt.ylabel("anomaly (ppm)")
plt.title("Seasonal Cycle of cfc113")
cbar.set_label("year")
plt.legend()
            
            
        


