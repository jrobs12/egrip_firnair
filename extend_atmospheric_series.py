# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 15:49:29 2025

@author: jrobs

Introduction
------------------------------------------
Final Workflow for correcting the Atmospheric Histories to use as a Firn Model Input. This notebook works to
extend the dataset used for the NEEM analysis in 2009. While the workflow is optimised for the EGRIP 2018 Greenland 
Project, the codecan be easily modified to be used for the Muller Project in 2025. The same workflow is used for all 
species, with some exceptions where data is missing. These exceptions are documented in the report as well as in comments
in the code itself. We found that for CFC113 there existed a gap between 2014 and 2016 in the monthly resolution data
which we filled with weekly resolution data re-sampled to monthly.

Gathered Data
-------------------------------------------

Data is gathered from the NOAA database from the available arctic stations in monthly resolution; alert (ALT), barrow 
(BRW) and summmit (SUM). Alert and Barrow are corrected to become an inferred Summit dataset, since Summit is located 
nearest the EGRIP sampling site. We used both flask and insitu data since it showed to be very similar. 

Workflow Breakdown
--------------------------------------------

(1) Data is read in and if applicable, in situ and flask data are merged. We start the data at 2000 to allow for adequete
    overlap between the new data and NEEM.
    
(2) We take a 12 month rolling mean to removed the seasonal cycle of the data.

(3) We subtract the rolling mean and find the mean seasonal cycle at each station in order to add this back in later
    steps
    
(4) Calculate the offset between either Barrow and Alert and the reference station (Summit), this is done by looking
    at the mean difference between the smoothed data (eg. rolling mean at sum - rolling mean at alt). If the found 
    offset is no better than 0 we take an offset of 0.

(5) Calculate the ratio between the stations and reference by minimising the least squares regression and finding the 
    geometric slope between the datasets. If this slope is not significantly different than 1 at the 0.01 level, we take 
    a slope of 1.
    
(6) Apply the coefficients to the station data and merge with the NEEM dataset. Calculate the uncertainty in our new
    measurements as the mean standard deviation between the residuals of the new corrected datasets.
    
    
Figures
-------------------------------------------
The code will produce 7-8 figures for each species. If there exists both flask and insitu data for a particular species
the code will plot both types of data to ensure that they are sufficiently similar to be mixed together. 

figure 1:
    Raw data at each station with the 12 month rolling mean.

figure 2:
    mean seasonal cycle at each station with the actual values scattered. We look for if there are big differences in
    seasonal cycle between stations at different parts of the year.

figure 3:
    The residuals between the reference station and the other stations. Both the residuals between the rolling means
    and the residuals between the de-trended (straight) data are plotted. We are looking for a resiudal with no trend
    ideally centred around 0 .
    
figure 4:
    Plot the smooth data against the reference station. Line of best fit is printed on the plot and both regression
    lines are plotted (a v b and b v a)
    
figure 5:
    The corrected data plotted with the NEEM dataset and Summit. Assess the overlap between the datasets and make sure
    that they look the same.

figure 6:
    final dataset from 2000 to 2019, check that the uncertainties are okay and that the dataset has no clear jump between
    the NEEM time and the new data.


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
import scipy.interpolate as interp
import scipy.linalg as alg #for linear regression
import datetime
import scipy.stats as stats
import warnings

#switch working directory to the USRA_2025, this is specifically because my computer hates me
cwd = os.getcwd()
wd = 'C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Code'
final_data_directory = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\Final\\"
if cwd != wd:
    os.chdir(wd)

#which columns to grab for each type of data
#want the year, month and value (+std if available)
gg_columns_flask = (1,2,3) #the first column is the station_id and throws an error when read-in
gg_columns_insitu = (1,2,10,11) #this data is in a different format, grab the same columns + the standard deviations
hal_columns = (0,1,2,3,4) # for the halocompounds since those are also in a different form

#the only warnings that I'm getting are from pandas and I can't get rid of them no matter what I try
warnings.filterwarnings("ignore", category=UserWarning)

#%% define calculation functions

#add date column
def add_dates(df):
    """Function to add a decimal date column to a set of data, given that the data has the year and the month.
    
       The function will calculate the decimal date as year + (dayofyear-1)/365, and assumes that the dayofyear 
       the 15th of each month
    
    Inputs
    --------------------------
    df: pandas dataframe
        dataframe containing the year and month column
    returns
    --------------------------
    df: pandas dataframe
        dataframe with one extra column containing the decimal dates"""
    df["date"] = 0
    df["day"] = 15 # add a day column
    df["dt"] = pd.to_datetime(df[["year", "month", "day"]])
    df["date"] = df["year"].astype("float") + df["month"].astype("float")/12
    df["date"] = df["dt"].dt.year + (df["dt"].dt.dayofyear-1)/365
    df.drop(labels = ["day", "dt"], axis = 1, inplace = True) #clean up the data, axis = 1 is columnwise

    return df

#tidy the data by cropping to relevant dates
def clean_data(df, start_time = 2000, end_time = 2019):
    """Add date column using the add_dates function, standardize dates to start at 2000 and end in 2019, set index of 
    dataframe to date. 
    inputs
    -----------
    df: pandas dataframe
        raw dataframe that was read in
    start_time: float, default 2000
        the year at which you want to start every timeseries
    end_time: float, default 2019
        The year at which you want to end the timeseries
    
    returns
    -------------
    df: pandas dataframe
        the cleaned dataframe with date column and new index
    """
    df = add_dates(df) #get dates in clean format
    #df = find_seasonal_cycle(df, plot = False) #get anomalies
    df = df[df["date"] > start_time] #grab all of the data that NEEM doesn't have with 10 years overlap
    df = df[df["date"] < end_time] #don't need all the way to 2023 for EGRIP project
    df.set_index(["date"], inplace = True, drop = False)
    
    return df

#resample weekly to monthly data
def resample_monthly(raw_data, species):
   """"resample weekly data to monthly at 3 arctic sites
   
   inputs
   -----------------
   raw_data: pandas dataframe
       data downloaded as discrete data from NOAA database, contains all of the weekly data at all stations
   species: string
       the name of the chemical species as listed in the raw data
      
        
    returns
    ------------------
    alt_flask: pandas dataframe
        resampled data at monthly resolution at the alert station
    brw_flask: pandas dataframe
        reseampled data at monthly resolution at the barrow station
    sum_flask: pandas dataframe
        resampled data at monthly resolution at the summit station
   """
    
   
   raw_alt_flask = raw_data[raw_data["site"] == "alt"]
   raw_brw_flask = raw_data[raw_data["site"] == "brw"]
   raw_sum_flask = raw_data[raw_data["site"] == "sum"]
  
    #find the monthly means of the data
   sample_data = [raw_alt_flask, raw_brw_flask, raw_sum_flask]
   for s, d in enumerate(sample_data):
       d.loc[:,"value"] = d[species]
       d.set_index(np.arange(len(d["value"])), inplace = True)
       d["year"] = 0
       d["month"] = 0
       d.loc[:, "date"] = pd.to_datetime(d["yyyymmdd hhmmss"].astype("datetime64[ns]"))
       for i in range(len(d["value"])): #add year and month column
           d.loc[i, 'year'] = d["date"][i].year
           d.loc[i, "month"] = d["date"][i].month
       final_data = d[d["year"] >= 2000][["year", "month", "value", "date"]]
       final_data_mean = final_data.groupby(["year", "month"], as_index = False).mean()
 
       
       if s == 0:
           alt_flask = final_data_mean
       
       if s == 1:
           brw_flask = final_data_mean
       
       if s == 2:
           sum_flask = final_data_mean
           
   return alt_flask, brw_flask, sum_flask

#get all of the data into one big dataframe
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
    #df.set_index(["date"], drop = False, inplace = True)
    df["year"] = year
    df["month"] = month
    df["alt"] = alt_data["value"]
    df["brw"] = brw_data["value"]
    df["sum"] = sum_data["value"]
    
    
    return df

def calculate_residuals(df, stations = ["alt", "brw", "sum"], reference_station = "sum"):
    """Finds residuals for smoothed data (rolling mean) and data with the seasonal cycle removed. 
    
    (1) The function will first calculate the 12 month rolling mean over the data, enforcing 12 months of non-nan 
    data to produce a valid point.
    
    (2) It then calculates the anomaly in reference to this rolling mean by subtracting the rolling mean from the raw data.
    
    (3) Finds the mean seasonal cycle at each station.
    
    (4) Subtracts the mean seasonal cycle from the raw data.
    
    (5) Finds the residual of the data without the seasonal cycle
   ---------------
   df: pandas dataframe
       dataframe produced using the make_data_file function. A N x 6 dataframe with a date, year and month column
       along with monthly average data for alert, barrow and summit stations in the last 3 columns.
       
    stations: list, default: ['alt', 'brw', 'sum']
        the stations that the data was collected from. Also the names of the data columns in df
        
    reference_station: string, default: "sum"
        The station that the residuals will be calculated in reference to
        
    Returns
    --------------
    The inputted dataframe with 3 extra columns for each station
    columns ending in:
        _rolling    -> the 12 month rolling mean of the data
        _anom       -> the raw data - 12 month rolling mean
        _res        -> summit_rollingmean - station_rollingmean
        _straight    -> the raw data minus the mean seasonal cycle
        _straight_res -> the residual of the data without the seasonal cycle
    """
    
    for s in stations:
        df[s+"_rolling"] = df[s].rolling(12, min_periods = 12, center = True).mean() #find the rolling mean
        df[s+"_anom"] = df[s] - df[s+"_rolling"] #calculate the anomaly
        df_seasonal_means = df.groupby("month").mean()
        
        # station_data = df[[s, "year", "month"]]
        # df[s + "_anom"] = 0 #make anomaly column
        
        # df_yearly = station_data.groupby('year').mean() #find the yearly mean for that year and station
    
        # for y in df_yearly.index: #year
        #     value_at_year = station_data[station_data["year"]==y][s]
        #     if (value_at_year/value_at_year).sum() == 12: 
        #         #check to see if we have a full year. Nans are not counted as numbers so the test will fail
        #         df.loc[df["year"]==y, s+"_anom"] = value_at_year-df_yearly[s][y] 
        #         #subtract to make anomaly
        #     else: df.loc[df["year"]==y, s+"_anom"] = np.nan #make the year nans if the dataset is not complete
    
    
    for s in stations:
        df[s+"_rolling_res"] = df[reference_station+"_rolling"] - df[s+"_rolling"]
        #find the residuals on the rolling means (for offset)
        
        #subtract the mean monthly cycle
        station_data = df[["year", "month", s]]
        df[s+"_straight"] = 0
        for m in df_seasonal_means.index: #months
            cycle_at_month = df_seasonal_means[s+"_anom"][m]
            data_at_month = station_data[station_data["month"] == m]
            df.loc[df["month"]==m, s+"_straight"] = data_at_month[s] - cycle_at_month
    for s in stations:
        df[s+'_straight_res'] = df[reference_station+"_straight"] - df[s +"_straight"] 
        #calculate the residual between the cleaner data

    return df

#find the offset and the ratio to apply to data
def find_coefficients(df, stations = ["alt", "brw"], alt = True):
    """Calculates the Offset between 2 stations, and finds the correct ratio to apply in order to correct the data.
    (1) Determines if the mean of the difference between the reference station and the other stations is close enough to 0
        to enforce an offset of 0. 
            if the mean of the residuals is less than the standard error we say that there is a 0 offset.
            
    (2) Calculates the geometric slope between the data without the trend line in order to find the ratio to apply.
            Solves the least squares optimisation problem ax = b where a is one set of data and b is the other.
            
    (3) determines if the slope is significantly different than one using the residues of the optimisation.
        if the slope is significant at a 1% level
        
    inputs
    --------------
    df: pandas dataframe
        dataframe on which to compute the coefficients, must contain the residual of the rolling mean, and of the detrended
        data
    stations: list of strings, default ["alt", "brw"] 
        the station on which you want to calculate the coefficients
    alt: Boolean, default True
        select true if there is data from an alert station to be analysed
        
    returns
    -----------------
    all_coefficients: dictionary
        dictionary containing a dictionary for each station that contains the coefficients """
    all_coefficients = {}
    
    if alt == False: #for when we don't have data from alert
        df.drop(labels = "alt", inplace = True, axis = 1)
        
    for s in stations:
        station_coefficients = {} #make empty dictionary for each station

        #ideally the mean of data without an offset is 0, determine if the mean is close enough to 0 to ignore it
        mean_res_rolling = df[s + "_rolling_res"].mean()
        
        n_rolling = len(df[s + "_rolling_res"])
        
        std_rolling = df[s + "_rolling_res"].std()
        
        #find the standard error
        se = std_rolling / np.sqrt(n_rolling)
        print(se, mean_res_rolling)
        
        #check if different than 0
        if np.abs(mean_res_rolling) < (std_rolling / np.sqrt(n_rolling)):
            station_coefficients["offset"] = 0
        else:
            offset = -1*mean_res_rolling #find offset between data
            print(offset)
            station_coefficients["offset"] = offset
        
        #find the geometric slope of the data using the residuals without the trend
        
        df = df.dropna() #need to drop the nans to run scipy.linalg.lstsq
        
        #find correlation between the stations
        r, p = stats.pearsonr(df["sum_anom"], df[s+"_anom"])
        print(f"{s} Correlation: r = {r}, p = {p}")
        
        #reshape the data to be an Nx1 array
        sum_reshape = df["sum_anom"].values.reshape(df["sum_anom"].shape[0], -1)
        station_reshape = df[s+"_anom"].values.reshape(df[s+"_anom"].shape[0], -1)
        
        sum_v_station_rolling = alg.lstsq(sum_reshape, df[s+"_anom"])[0] #find the slope
        station_v_sum_rolling = alg.lstsq(station_reshape, df["sum_anom"])[0] #find the slope
        
        station_coefficients["sum_v_station"] = sum_v_station_rolling
        station_coefficients["station_v_sum"] = station_v_sum_rolling
        
        #we want the geometric slope such that x is summit and y is the station in ax = b

        slope = 0.5*(sum_v_station_rolling + 1/station_v_sum_rolling)
        residual_corrected = (df["sum_anom"] - 1/slope * df[s+"_anom"]).std()
        residual = (df["sum_anom"] - df[s+"_anom"]).std()
        print(f"Residual After Correction s: {residual_corrected}, vs {residual}")
        
        if p < 0.01:
            #we want the geometric slope such that x is summit and y is the station in ax = b
            station_coefficients["ratio"] = slope #append to dictionary
        else:
           station_coefficients["ratio"] = 1
            
        all_coefficients[s] = station_coefficients #append each set of coefficients to the big dictionary
    return all_coefficients

def do_correction(df, coefficients, stations = ["alt", "brw"]):
    """apply the correction coefficients to find an inferred summit timeseries. The function will print the
    ratio and offset that will be applied to the data along with the station.
    
    inputs
    -----------------
    df: pandas dataframe
        the dataframe that contains the data that will be corrected
    coefficients: dictionary
        a nested dictionary of stations and their coefficients
    stations: list of strings, default is ["alt", "brw"] 
        list of stations that are to be corrected
    
    returns
    --------------------
    df: pandas dataframe
        the input dataframe with 6 new columns containing
        
        s + _inf_rolling -> the rolling mean with the offset applied
        s + _inf_rolling_res -> the residual of the s+ _inf_rolling with summit data
        s + _inf -> the rolling mean with the offset applied as well as the detrended data with the ratio applied
        """
    
    for s in stations:
        ratio = coefficients[s]["ratio"]
        offset = coefficients[s]["offset"]
        
        df[s+"_inf_rolling"] = (df[s+"_rolling"] - offset) #find the new rolling mean when the offset is added
        df[s+"_inf_rolling_res"] = df["sum_rolling"] - df[s+"_inf_rolling"] #find the residual
        df[s+"_inf"] = (df[s+"_rolling"] - offset) + df[s+"_anom"]/ratio #add wiggles back onto the data
        df[s+"_inf_res"] = df["sum"] - df[s+"_inf"] #find the new residual of the data
        
        print(f" At {s} the offset applied is {offset} and the ratio applied  = {1/ratio}")
        
    return df

def make_final_dataset(df, df_NEEM, final_path, species, alt = True):
    """ take the average of all the data, fill in any gaps using the summit seasonal cycle and find the standard error
    
    inputs
    -------------------
    df: pandas dataframe
        dataframe containing the corrected data
    df_NEEM: pandas dataframe
        dataframe containing the columns date, value, std that were used in the NEEM analysis
    final_path: string
        the location on the computer that the final dataset will be saved to
    species: string
        the name of the chemical species IN ALL CAPS that will be added to the name of the file
    alt: boolean, default = True
        select true if there exists data from the alert station
        
    returns
    --------------------
    df: pandas dataframe
        dataframe that contains the uncertainty and merged value for the data at EGRIP from start_time to end_time
    egrip_dataset: pandas dataframe
        dataframe containing the atmospheric history from pre-start_time to end_time"""
    if alt == True:
        final_data_labels = ["alt_inf", "brw_inf", "sum", 'NEEM']
        final_std_labels = ["alt_inf_rolling_res", "brw_inf_rolling_res"]
    else:
        final_data_labels = ["brw_inf", "sum", "NEEM"]
        final_std_labels = ["brw_inf_rolling_res"]
    
    #calculate uncertainty on the new datat0o9iiuyg
    uncertainty_max = np.nanmax(df[final_std_labels].values)
    uncertainty_min = np.nanmin(df[final_std_labels].values)
    uncertainty = uncertainty_max - uncertainty_min #take the difference between the two inferred datasets
    
    df["std"] = uncertainty #average difference between the residuals
    #take std between the max and min differences to get one number
    
    df_NEEM.set_index(["date"], drop = False, inplace = True)
    end_NEEM = df_NEEM["date"].values[-1]
    NEEM_overlap = df_NEEM[df_NEEM["date"] >= df["date"].values[0]] #find where our data starts in the NEEM dataset
    df_no_overlap = df[df["date"] > end_NEEM]
    extra_nans = pd.Series(np.ones(len(df_no_overlap.index))*np.nan, index = df_no_overlap.index)
    #use these to pad out the NEEM for merging
    
    df["NEEM"] = pd.concat([NEEM_overlap["value"], extra_nans])
    df["NEEM_unc"] = pd.concat([NEEM_overlap["std"], extra_nans])
    
    df["value"] = df[final_data_labels].mean(axis = 1)
    df["std"] = df[["std", "NEEM_unc"]].mean(axis = 1) #average the NEEM uncertainties with ours
    
    NEEM_no_overlap = df_NEEM[df_NEEM["date"] <= df["date"].values[0] + 196/365] #because the rolling mean centres
    #the window, we need to add 6 months since we don't have inferred data for those times
    
    egrip_dataset = pd.concat([NEEM_no_overlap, df[["date", "value", "std"]]], ignore_index= True)
    
    #save the file
    egrip_dataset.to_csv(final_path+f"SCENARIO_EGRIP18_{species}.csv")
    return df, egrip_dataset

#%% define the plotting functions

#define colormap
cmap = mpl.colormaps["berlin"]
cmap2 = sns.color_palette(palette='icefire')
n_lines = 7
colors = cmap(np.linspace(0,1,n_lines))

#flask and insitu data check
def plot_figure0(df_flask, df_insitu, fig_num, species, station):
    """plot the overlap of flask and insitu data at stations where both exist"""
        

    plt.figure(fig_num, clear = True)
    plt.plot(df_flask["date"], df_flask["value"], color = colors[3], label = "Flask Data")
    plt.plot(df_insitu["date"], df_insitu["value"], color = colors[6], label = "Insitu Data")
    plt.xlabel("time"), plt.ylabel("concentration (ppm)")
    plt.title(f"{species} data from 2 sources at {station}")
    plt.legend()

#raw data
def plot_figure1(df, fig_num, species, final_path, alt = True):
    """plot the raw data from the station columns"""
    plt.figure(fig_num, clear = True) 
    if alt == True:
        plt.plot(df["date"], df["alt"], label = "ALERT", color = cmap2[0], ls = "dashed")
    plt.plot(df["date"], df["brw"], label = "BARROW", color = cmap2[4], ls = "dashed")
    plt.plot(df["date"], df["sum"], label = "SUMMIT", color = cmap2[5], ls = "dashed")
    
    if alt == True:
        plt.plot(df["date"], df["alt_rolling"], color = cmap2[0])
    plt.plot(df["date"], df["brw_rolling"], color = cmap2[4])
    plt.plot(df["date"], df["sum_rolling"], color = cmap2[5])

    plt.xlabel("time")
    plt.ylabel(f"{species} concentration (ppm)")
    plt.title(f"Atmospheric Concentration of {species} in the High Arctic")
    plt.legend(title = "--Station--")

    plt.savefig(final_path + "figure1.png")

#seasonal cycle
def plot_figure2(df, fig_num, species, final_path, alt = True):
    """plot the seasonal cycle of a species at different stations"""
    plt.figure(fig_num, clear = True)
    #find seasonal means at each station
    df_seasonal_means = df.groupby("month").mean()
    
    #plot the mean seasonal cycle
    if alt == True:
        plt.plot(df_seasonal_means.index, df_seasonal_means["alt_anom"], color = cmap2[0], label = "alert")
    plt.plot(df_seasonal_means.index, df_seasonal_means["brw_anom"], color = cmap2[4], label = "barrow")
    plt.plot(df_seasonal_means.index, df_seasonal_means["sum_anom"], color = cmap2[5], label = "summit")
    
    #plot each individual point
    if alt == True:
        plt.scatter(df["month"], df["alt_anom"], color = cmap2[0], alpha = 0.5)
    plt.scatter(df["month"], df["brw_anom"], color = colors[4], alpha = 0.5)
    plt.scatter(df["month"], df["sum_anom"], color = colors[5], alpha = 0.5)

    plt.xlabel("month")
    plt.ylabel("seasonal change from mean (ppm)")
    plt.title(f"Seasonal Cycle of {species} at 3 stations in the High Arctic")
    plt.legend(loc = "lower left", title = "--Station--")

    plt.savefig(final_path + "figure2.png")
    
#residuals
def plot_figure3(df, fig_num, species, final_path, alt = True):
    """plot the residuals of the detrended data (_anom) and the rolling means (_rolling)"""
    if alt == True:
        f = plt.figure(fig_num, clear = True)
        axs = f.subplots(1,2)
        
        ax = axs[0]
        ax.plot(df["date"], df["alt_straight_res"], color = cmap2[0], ls = "dashed")
        ax.plot(df["date"], df["alt_rolling_res"], color = cmap2[0], label = "alert")
        
        ax.plot(df["date"], df["sum_straight_res"], color = cmap2[5], ls = "dashed")
        ax.plot(df["date"], df["sum_rolling_res"], color = cmap2[5], label = "summit")
        
        ax.legend(title = "-- Station --")
        ax.set_xlabel("time")
        ax.set_ylabel("residual (station - summit)")
        ax.set_title(f"Residuals of {species} ast Alert \n -- Detrended Resdiuals   -Rolling Mean Residuals")


        ax = axs[1]
        ax.plot(df["date"], df["brw_straight_res"], color = cmap2[4], ls = "dashed")
        ax.plot(df["date"], df["sum_straight_res"], color = cmap2[5], ls = "dashed")
        
        ax.plot(df["date"], df["brw_rolling_res"], color = cmap2[4], label = "barrow")
        ax.plot(df["date"], df["sum_rolling_res"], color = cmap2[5], label = "summit")
        
        ax.legend(title = "-- Station --")
        ax.set_xlabel("time")
        ax.set_ylabel("residual (station - summit)")
        ax.set_title(f"Residuals of {species} at Barrow \n -- Detrended Resdiuals   -Rolling Mean Residuals")
    
    else:
        f = plt.figure(fig_num, clear = True)
        plt.plot(df["date"], df["brw_straight_res"], color = cmap2[4], ls = "dashed")
        plt.plot(df["date"], df["sum_straight_res"], color = cmap2[5], ls = "dashed")
        
        plt.plot(df["date"], df["brw_rolling_res"], color = cmap2[4], label = "barrow")
        plt.plot(df["date"], df["sum_rolling_res"], color = cmap2[5], label = "summit")
        
        plt.legend(title = "-- Station --")
        plt.xlabel("time")
        plt.ylabel("residual (station - summit)")
        plt.title(f"Residuals of {species} at Barrow \n -- Detrended Resdiuals   -Rolling Mean Residuals")
    
#linear regression    
def plot_figure4(df, coefficients, fig_num, species, final_path, alt = True):
    if alt == True:
        f = plt.figure(fig_num, clear = True)
        axs = f.subplots(1,2)
        
        alt_slope = coefficients["alt"]["ratio"]
        brw_slope = coefficients["brw"]["ratio"]
        
        ax = axs[0]
        ax.plot(df["sum_anom"], alt_slope*df["sum_anom"], color = cmap2[0], label = "alert line")
        ax.plot(df["sum_anom"], coefficients["alt"]["sum_v_station"]*df["sum_anom"], color = cmap2[0], ls = "dotted" )
        ax.plot(df["sum_anom"], 1/coefficients["alt"]["station_v_sum"]*df["sum_anom"], ls = "dotted")
        
        ax.scatter(df["sum_anom"], df["alt_anom"], color = cmap2[0], alpha = 0.7)
        
        ax.set_xlabel("summit")
        ax.set_ylabel("alert")
        ax.set_title("Regression of  detrended Summit Against Alert")
        ax.legend()
        ax.grid()
        
        ax = axs[1]
        ax.plot(df["sum_anom"], brw_slope*df["sum_anom"], color = colors[4], label = "barrow line")
        ax.plot(df["sum_anom"], coefficients["brw"]["sum_v_station"]*df["sum_anom"], ls = "dashed", color = colors[4])
        ax.plot(df["sum_anom"], 1/coefficients["brw"]["station_v_sum"]*df["sum_anom"], ls = "dashed", color = colors[4])
    
    
        ax.scatter(df["sum_anom"], df["brw_anom"], color = cmap2[4], alpha = 0.7)
        
        ax.set_xlabel("summit")
        ax.set_ylabel("barrow")
        ax.legend()
        ax.set_title("Regression of detrended Summit Against Barrow")
        ax.grid()
    
    else:
        plt.figure(fig_num, clear = True)
        brw_slope = coefficients["brw"]["ratio"]
        plt.plot(df["sum_anom"], brw_slope*df["sum_anom"], color = colors[4], label = "barrow line")
        plt.plot(df["sum_anom"], coefficients["brw"]["sum_v_station"]*df["sum_anom"], ls = "dashed", color = colors[4])
        plt.plot(df["sum_anom"], 1/coefficients["brw"]["station_v_sum"]*df["sum_anom"], ls = "dashed", color = colors[4])
    
    
        plt.scatter(df["sum_anom"], df["brw_anom"], color = cmap2[4], alpha = 0.7)
        
        plt.xlabel("summit")
        plt.ylabel("barrow")
        plt.legend()
        plt.title("Regression of detrended Summit Against Barrow")
        plt.grid()
    
#corrected data
def plot_figure5(df, coefficients, fig_num, species, final_path, alt = True):
    if alt == True:
        f = plt.figure(fig_num, clear = True)
        axs = f.subplots(2,2)
        
        #plot alert data
        ax = axs[0,0]
        ax.plot(df["date"], df["alt_inf"], color = cmap2[0], label = "Inferred Summit")
        ax.plot(df["date"], df["sum"], color = cmap2[5], label = "Actual Summit")
        
        ax.plot(df["date"], df["alt"], ls = "dashed", color = cmap2[0])
        ax.plot(df["date"], df["alt_rolling"], color = cmap2[0], ls = "-.")
        ax.plot(df["date"], df["sum_rolling"], color = cmap2[5])
        ax.plot(df["date"], df["alt_rolling"] - coefficients["alt"]["offset"], color = cmap2[0])
        
        ax.set_xlabel("time")
        ax.set_ylabel(f"{species} concentration (ppm)")
        ax.set_title(f"Effects of Correction on {species} Data at ALERT")
        ax.legend()
        
        #plot barrow data
        ax = axs[0,1]
        ax.plot(df["date"], df["brw_inf"], color = cmap2[4], label = "Inferred Summit")
        ax.plot(df["date"], df["sum"], color = cmap2[5], label = "Actual Summit")
        
        ax.plot(df["date"], df["brw"], ls = "dashed", color = cmap2[4])
        ax.plot(df["date"], df["brw_rolling"], color = cmap2[4], ls = "-.")
        ax.plot(df["date"], df["sum_rolling"], color = cmap2[5])
        ax.plot(df["date"], df["brw_rolling"] - coefficients["brw"]["offset"], color = cmap2[4])
        
        ax.set_xlabel("time")
        ax.set_ylabel(f"{species} concentration (ppm)")
        ax.set_title(f"Effects of Correction on {species} Data at BARROW")
        ax.legend()
        
        #plot alert residuals
        ax = axs[1,0]
        ax.plot(df["date"], df["alt_inf_res"], color = cmap2[0], ls = "dashed", label = "Residuals")
        ax.plot(df["date"], df["alt_inf_rolling_res"], color = cmap2[0], label = "Residual on Rolling Mean")
        ax.set_xlabel("time")
        ax.set_ylabel(f"(sum - alt) (ppm)")
        ax.set_title(f"Residuals after correction at ALERT")
        ax.axhline(color = cmap2[5], ls = "dotted")
        ax.grid()
        ax.legend()
        
        ax = axs[1,1]
        ax.plot(df["date"], df["brw_inf_res"], color = cmap2[4], ls = "dashed", label = "Residuals")
        ax.plot(df["date"], df["brw_inf_rolling_res"], color = cmap2[4], label = "Residual on Rolling Mean")
        ax.set_xlabel("time")
        ax.set_ylabel(f"(sum - brw) (ppm)")
        ax.axhline(color = cmap2[5], ls = "dotted")
        ax.set_title(f"Residuals after correction at BARROW")
        ax.grid()
        ax.legend()
        
    else:
        f = plt.figure(fig_num, clear = True)
        axs = f.subplots(2,1)
        
        #plot just the barrow data
        ax = axs[0]
        ax.plot(df["date"], df["brw_inf"], color = cmap2[4], label = "Inferred Summit")
        ax.plot(df["date"], df["sum"], color = cmap2[5], label = "Actual Summit")
        
        ax.plot(df["date"], df["brw"], ls = "dashed", color = cmap2[4], label = "Uncorrected Barrow")
        ax.plot(df["date"], df["brw_rolling"], label = "Barrow Rolling Mean", color = cmap2[4])
        ax.plot(df["date"], df["sum_rolling"], label = "Summit Rolling Mean", color = cmap2[5])
        ax.plot(df["date"], df["brw_rolling"] - coefficients["brw"]["offset"], 
                label = "Barrow minus offset", color = cmap2[3])
        
        ax.set_xlabel("time")
        ax.set_ylabel(f"{species} concentration (ppm)")
        ax.set_title(f"Effects of Correction on {species} Data at BARROW")
        ax.legend()
        
        #plot the residuals
        ax = axs[1,0]
        ax.plot(df["date"], df["brw_inf_res"], color = cmap2[4], ls = "dashed", label = "Residuals")
        ax.plot(df["date"], df["brw_inf_rolling_res"], color = cmap2[4], label = "Residual on Rolling Mean")
        ax.set_xlabel("time")
        ax.set_ylabel(f"(sum - brw) (ppm)")
        ax.axhline(color = cmap2[5], ls = "dotted")
        ax.set_title(f"Residuals after correction at BARROW")
        ax.grid()
        ax.legend()
        

#NEEM overlap
def plot_figure6(df, df_NEEM, fig_num, species, final_path, alt = True):
    """Plot all available data between 2000 and 2019"""
    NEEM_period = df_NEEM[df_NEEM["date"] >= 2000]
    
    plt.figure(fig_num, clear = True)
    
    plt.plot(NEEM_period["date"], NEEM_period["value"], color = cmap2[2], label = "NEEM")
    if alt == True:
        plt.plot(df["date"], df["alt_inf"], color = cmap2[0], label = "Alert")
    plt.plot(df["date"],df["brw_inf"], color = cmap2[4], label = "Barrow")
    plt.plot(df["date"],df["sum"], color = cmap2[5], label = "Summit")
    
    plt.xlabel("time")
    plt.ylabel("concentration (ppm)")
    plt.title(f"Concentration of {species} in the High Arctic")
    
def plot_figure7(df, df_NEEM, fig_num, species, final_path):
    """plot the final dataset with errorbars"""
    last_NEEM = df_NEEM["date"].values[-1]
    plt.figure(fig_num, clear = True)
    plt.plot(df["date"], df["value"], color = cmap2[3])
    plt.fill_between(df["date"], df["value"] + df["std"], df["value"] - df["std"], color = cmap2[3], alpha = 0.4)
    plt.axvline(last_NEEM, color = colors[5])
    
    plt.xlim(1990)
    plt.xlabel("time")
    plt.ylabel("concentration (ppm)")
    plt.title(f"Concentration of {species} in the High Arctic \n Merged with NEEM")
#%% test with co2
co2_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\co2\\*.txt" #where the data is stored
co2_files = glob.glob(co2_path)
co2_path = final_data_directory+"co2\\" #where to save the figures when we're done

#load in each station's data
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

#average barrow together
co2_brw = co2_brw_flask
co2_brw["data_insitu"] = co2_brw_insitu["value"]

#check that both barrows look similar
plot_figure0(co2_brw_flask, co2_brw_insitu, 0, "co2", "BARROW")


#make datafile and find all of the values
co2_data = make_data_file(co2_alt_flask, co2_brw, co2_sum_flask, co2_brw["date"], co2_brw["month"], co2_brw["year"])
#find the seasonal cycle, subtract it and look at the residuals
co2_data = calculate_residuals(co2_data)

plot_figure1(co2_data, 1, "co2", co2_path) #plot the raw values and the rolling mean
plot_figure2(co2_data, 2, "co2", co2_path) #plot the seasonal cycle at each station

co2_coefficients = find_coefficients(co2_data) #find coefficients to apply to make the data
plot_figure3(co2_data, 3, "co2", co2_path) #plot the residuals (rolling and detrended)
plot_figure4(co2_data, co2_coefficients, 4, "co2", co2_path) #check that the best fit line actually has 0 intercept
#and is fitting the data well
co2_data = do_correction(co2_data, co2_coefficients) #apply the offset and ratio
plot_figure5(co2_data, co2_coefficients, 5, "co2", co2_path) #check that the correction worked
plot_figure6(co2_data, co2_NEEM, 6, "co2", co2_path) #check the overlap with NEEM

#average the dataset together for the final data
co2_data, co2_egrip = make_final_dataset(co2_data, co2_NEEM, final_data_directory, "CO2")
#plot the final dataset
plot_figure7(co2_egrip, co2_NEEM, 7, "co2", co2_path)

#%% CH4 Data
ch4_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\ch4\\*.txt"
ch4_files = glob.glob(ch4_path)
ch4_path = final_data_directory + "ch4\\"

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


ch4_brw = ch4_brw_flask
ch4_brw["data_insitu"] = ch4_brw_insitu["value"]

#check that both barrows look similar
plot_figure0(ch4_brw_flask, ch4_brw_insitu, 10, "ch4", "BARROW")

#make datafile and find all of the values
ch4_data = make_data_file(ch4_alt_flask, ch4_brw, ch4_sum_flask, ch4_brw["date"], ch4_brw["month"], ch4_brw["year"])
#find the seasonal cycle, subtract it and look at the residuals
ch4_data = calculate_residuals(ch4_data)

plot_figure1(ch4_data, 11, "ch4", ch4_path) #plot the raw values and the rolling mean
plot_figure2(ch4_data, 12, "ch4", ch4_path) #plot the seasonal cycle at each station

ch4_coefficients = find_coefficients(ch4_data) #find coefficients to apply to make the data
plot_figure3(ch4_data, 13, "co2", ch4_path) #plot the residuals (rolling and detrended)
plot_figure4(ch4_data, ch4_coefficients, 14, "ch4", ch4_path) #check that the best fit line actually has 0 intercept
#and is fitting the data well
ch4_data = do_correction(ch4_data, ch4_coefficients) #apply the offset and ratio
plot_figure5(ch4_data, ch4_coefficients, 15, "ch4", ch4_path) #check that the correction worked
plot_figure6(ch4_data, ch4_NEEM, 16, "ch4", ch4_path) #check the overlap with NEEM

#average the dataset together for the final data
ch4_data, ch4_egrip = make_final_dataset(ch4_data, ch4_NEEM, final_data_directory, "CH4")
#plot the final dataset
plot_figure7(ch4_egrip, ch4_NEEM, 17, "ch4", ch4_path)

#%% SF6 Data
sf6_path= "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\sf6\\*.txt"
sf6_files = glob.glob(sf6_path)
sf6_path = final_data_directory+"sf6\\"

sf6_alt_flask = pd.DataFrame( np.loadtxt(sf6_files[2], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
sf6_alt_flask = clean_data(sf6_alt_flask)


sf6_brw_flask = pd.DataFrame( np.loadtxt(sf6_files[3], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
sf6_brw_flask = clean_data(sf6_brw_flask)


sf6_brw_insitu = pd.DataFrame( np.loadtxt(sf6_files[0], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
sf6_brw_insitu = clean_data(sf6_brw_insitu)

sf6_brw = sf6_brw_flask
sf6_brw["data_insitu"] = sf6_brw_insitu["value"]


sf6_sum_flask = pd.DataFrame( np.loadtxt(sf6_files[4], usecols = gg_columns_flask), columns = 
                             ["year", "month", "value"])
sf6_sum_flask = clean_data(sf6_sum_flask)


sf6_sum_insitu = pd.DataFrame( np.loadtxt(sf6_files[5], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
sf6_sum_insitu = clean_data(sf6_sum_insitu)


sf6_brw = sf6_brw_flask
sf6_brw["data_insitu"] = sf6_brw_insitu["value"]

sf6_sum = sf6_sum_flask
sf6_sum["data_insitu"] = sf6_sum_insitu["value"]
#check that both barrows look similar
plot_figure0(sf6_brw_flask, sf6_brw_insitu, 29, "sf6", "BARROW")
plot_figure0(sf6_sum_flask, sf6_sum_insitu, 28, "sf6", "SUMMIT")

#average the barrows together
sf6_brw["mean_brw"] = sf6_brw[["value", "data_insitu"]].mean(axis = 1)
sf6_brw = sf6_brw[["mean_brw", "date", "year", "month"]]
sf6_brw.rename(columns = {"mean_brw" : "value"}, inplace = True)

#average the summits together
sf6_sum["mean_brw"] = sf6_sum[["value", "data_insitu"]].mean(axis = 1)
sf6_sum = sf6_sum[["mean_brw", "date", "year", "month"]]
sf6_sum.rename(columns = {"mean_brw" : "value"}, inplace = True)

sf6_NEEM = pd.DataFrame( np.loadtxt(sf6_files[1]), columns = 
                             ["date","value", "std"])

#make datafile and find all of the values
sf6_data = make_data_file(sf6_alt_flask, sf6_brw, sf6_sum_flask, sf6_brw["date"], sf6_brw["month"], sf6_brw["year"])
#find the seasonal cycle, subtract it and look at the residuals
sf6_data = calculate_residuals(sf6_data)

plot_figure1(sf6_data, 21, "sf6", sf6_path) #plot the raw values and the rolling mean
plot_figure2(sf6_data, 22, "sf6", sf6_path) #plot the seasonal cycle at each station

sf6_coefficients = find_coefficients(sf6_data) #find coefficients to apply to make the data
plot_figure3(sf6_data, 23, "sf6", sf6_path) #plot the residuals (rolling and detrended)
plot_figure4(sf6_data, sf6_coefficients, 24, "sf6", sf6_path) #check that the best fit line actually has 0 intercept
#and is fitting the data well
sf6_data = do_correction(sf6_data, sf6_coefficients) #apply the offset and ratio
plot_figure5(sf6_data, sf6_coefficients, 25, "sf6", sf6_path) #check that the correction worked
plot_figure6(sf6_data, sf6_NEEM, 26, "sf6", sf6_path) #check the overlap with NEEM

#average the dataset together for the final data
sf6_data, sf6_egrip = make_final_dataset(sf6_data, sf6_NEEM, final_data_directory, "SF6")
#plot the final dataset
plot_figure7(sf6_egrip, sf6_NEEM, 27, "sf6", sf6_path)

#%% CCL4 Data
ccl4_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\ccl4\\*.txt"
ccl4_files= glob.glob(ccl4_path)
ccl4_path = final_data_directory + "ccl4\\"
ccl4_brw_insitu = pd.DataFrame( np.loadtxt(ccl4_files[0], usecols= hal_columns), 
                               columns = ["year", "month", "value", "unc", "std"])
ccl4_brw_insitu = clean_data(ccl4_brw_insitu)

ccl4_sum_insitu = pd.DataFrame( np.loadtxt(ccl4_files[2], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
ccl4_sum_insitu = clean_data(ccl4_sum_insitu)

no_alt = pd.DataFrame(ccl4_brw_insitu["value"]*np.nan)
no_alt["date"] = ccl4_brw_insitu["date"]

ccl4_NEEM = pd.DataFrame( np.loadtxt(ccl4_files[1], ), columns = 
                             ["date","value", "std"])

#make datafile and find all of the values
ccl4_data = make_data_file(no_alt, ccl4_brw_insitu, ccl4_sum_insitu, ccl4_brw_insitu["date"], 
                           ccl4_brw_insitu["month"], ccl4_brw_insitu["year"])
#find the seasonal cycle, subtract it and look at the residuals
ccl4_data = calculate_residuals(ccl4_data, stations = ["brw", "sum"])

plot_figure1(ccl4_data, 31, "ccl4", ccl4_path, alt = False) #plot the raw values and the rolling mean
plot_figure2(ccl4_data, 32, "ccl4", ccl4_path, alt = False) #plot the seasonal cycle at each station

ccl4_coefficients = find_coefficients(ccl4_data, stations = ["brw"], alt = False) #find coefficients to apply to make the data

plot_figure3(ccl4_data, 33, "ccl4", ccl4_path, alt = False) #plot the residuals (rolling and detrended)
plot_figure4(ccl4_data, ccl4_coefficients, 34, "ccl4", ccl4_path, alt = False) #check that the best fit line actually has 0 intercept
#and is fitting the data well
ccl4_data = do_correction(ccl4_data, ccl4_coefficients, stations = ["brw"]) #apply the offset and ratio
plot_figure5(ccl4_data, ccl4_coefficients, 35, "ccl4", ccl4_path, alt = False) #check that the correction worked
plot_figure6(ccl4_data, ccl4_NEEM, 36, "ccl4", ccl4_path, alt = False) #check the overlap with NEEM

#average the dataset together for the final data
ccl4_data, ccl4_egrip = make_final_dataset(ccl4_data, ccl4_NEEM, final_data_directory, "CCL4", alt = False)
#plot the final dataset
plot_figure7(ccl4_egrip, ccl4_NEEM, 37, "ccl4", ccl4_path)

#%% CH3CCL3 Data
ch3ccl3_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\ch3ccl3\\*.txt"
ch3ccl3_files = glob.glob(ch3ccl3_path)
ch3ccl3_path = final_data_directory + "ch3ccl3\\"

ch3ccl3_brw_insitu = pd.DataFrame( np.loadtxt(ch3ccl3_files[0], usecols= hal_columns), 
                               columns = ["year", "month", "value", "unc", "std"])
ch3ccl3_brw_insitu = clean_data(ch3ccl3_brw_insitu)

ch3ccl3_brw_insitu1 = pd.DataFrame( np.loadtxt(ch3ccl3_files[1], usecols= hal_columns), 
                               columns = ["year", "month", "value", "unc", "std"])
ch3ccl3_brw_insitu1 = clean_data(ch3ccl3_brw_insitu1)

ch3ccl3_brw_insitu = pd.concat([ch3ccl3_brw_insitu, ch3ccl3_brw_insitu1], 
                               ignore_index= True)
ch3ccl3_brw_insitu = clean_data(ch3ccl3_brw_insitu)

ch3ccl3_sum_insitu = pd.DataFrame( np.loadtxt(ch3ccl3_files[2], usecols = hal_columns), columns = 
                             ["year", "month", "value", "unc", "std"])
ch3ccl3_sum_insitu = clean_data(ch3ccl3_sum_insitu)


ch3ccl3_NEEM = pd.DataFrame( np.loadtxt(ch3ccl3_files[3], ), columns = 
                             ["date","value", "std"])

#make datafile and find all of the values
ch3ccl3_data = make_data_file(no_alt, ch3ccl3_brw_insitu, ch3ccl3_sum_insitu, ch3ccl3_brw_insitu["date"], 
                           ch3ccl3_brw_insitu["month"], ch3ccl3_brw_insitu["year"])
#find the seasonal cycle, subtract it and look at the residuals
ch3ccl3_data = calculate_residuals(ch3ccl3_data, stations = ["brw", "sum"])

plot_figure1(ch3ccl3_data, 41, "ch3ccl3", ch3ccl3_path, alt = False) #plot the raw values and the rolling mean
plot_figure2(ch3ccl3_data, 42, "ch3ccl3", ch3ccl3_path, alt = False) #plot the seasonal cycle at each station

ch3ccl3_coefficients = find_coefficients(ch3ccl3_data, stations = ["brw"], alt = False) #find coefficients to apply to make the data

plot_figure3(ch3ccl3_data, 43, "ch3ccl3", ch3ccl3_path, alt = False) #plot the residuals (rolling and detrended)
plot_figure4(ch3ccl3_data, ch3ccl3_coefficients, 44, "ch3ccl3", ch3ccl3_path, alt = False) #check that the best fit line actually has 0 intercept
#and is fitting the data well
ch3ccl3_data = do_correction(ch3ccl3_data, ch3ccl3_coefficients, stations = ["brw"]) #apply the offset and ratio
plot_figure5(ch3ccl3_data, ch3ccl3_coefficients, 45, "ch3ccl3", ch3ccl3_path, alt = False) #check that the correction worked
plot_figure6(ch3ccl3_data, ch3ccl3_NEEM, 46, "ch3ccl3", ch3ccl3_path, alt = False) #check the overlap with NEEM

#average the dataset together for the final data
ch3ccl3_data, ch3ccl3_egrip = make_final_dataset(ch3ccl3_data, ch3ccl3_NEEM, final_data_directory, "CH3CCl3", alt = False)
#plot the final dataset
plot_figure7(ch3ccl3_egrip, ch3ccl3_NEEM, 47, "ch3ccl3", ch3ccl3_path)

#%% CFC11 Data
cfc11_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\cfc11\\*.txt"
cfc11_files = glob.glob(cfc11_path)
cfc11_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\Final\\cfc11\\"
cfc11_brw_insitu = pd.DataFrame(np.loadtxt(cfc11_files[1], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc11_brw_insitu1 = pd.DataFrame(np.loadtxt(cfc11_files[2], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc11_brw_insitu = pd.concat([cfc11_brw_insitu, cfc11_brw_insitu1], ignore_index = True)

cfc11_sum_insitu = pd.DataFrame(np.loadtxt(cfc11_files[4], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc11_sum_insitu = clean_data(cfc11_sum_insitu)
cfc11_brw_insitu = clean_data(cfc11_brw_insitu)

cfc11_NEEM = pd.DataFrame(np.loadtxt(cfc11_files[5], usecols = (0,1,2)), columns = 
                                ["date", "value", "std"])
cfc11_NEEM["value"] = cfc11_NEEM["value"]/1.0081

#make datafile and find all of the values
cfc11_data = make_data_file(no_alt, cfc11_brw_insitu, cfc11_sum_insitu, cfc11_brw_insitu["date"], 
                           cfc11_brw_insitu["month"], cfc11_brw_insitu["year"])
#find the seasonal cycle, subtract it and look at the residuals
cfc11_data = calculate_residuals(cfc11_data, stations = ["brw", "sum"])

plot_figure1(cfc11_data, 51, "cfc11",cfc11_path, alt = False) #plot the raw values and the rolling mean
plot_figure2(cfc11_data, 52, "cfc11", cfc11_path, alt = False) #plot the seasonal cycle at each station

cfc11_coefficients = find_coefficients(cfc11_data, stations = ["brw"], alt = False) #find coefficients to apply to make the data

plot_figure3(cfc11_data, 53, "cfc11", cfc11_path, alt = False) #plot the residuals (rolling and detrended)
plot_figure4(cfc11_data, cfc11_coefficients, 54, "cfc11", cfc11_path, alt = False) #check that the best fit line actually has 0 intercept
#and is fitting the data well
cfc11_data = do_correction(cfc11_data, cfc11_coefficients, stations = ["brw"]) #apply the offset and ratio
plot_figure5(cfc11_data, cfc11_coefficients, 55, "cfc11", cfc11_path, alt = False) #check that the correction worked

plot_figure6(cfc11_data, cfc11_NEEM, 56, "cfc11", cfc11_path, alt = False) #check the overlap with NEEM

#average the dataset together for the final data
cfc11_data, cfc11_egrip = make_final_dataset(cfc11_data, cfc11_NEEM, final_data_directory, "CFC11", alt = False)
#plot the final dataset
plot_figure7(cfc11_egrip, cfc11_NEEM, 57, "cfc11", cfc11_path)

#%% CFC12
cfc12_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\cfc12\\*.txt"
cfc12_files = glob.glob(cfc12_path)
cfc12_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\Final\\cfc12\\"
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

#make datafile and find all of the values
cfc12_data = make_data_file(no_alt, cfc12_brw_insitu, cfc12_sum_insitu, cfc12_brw_insitu["date"], 
                           cfc12_brw_insitu["month"], cfc12_brw_insitu["year"])
#find the seasonal cycle, subtract it and look at the residuals
cfc12_data = calculate_residuals(cfc12_data, stations = ["brw", "sum"])

plot_figure1(cfc12_data, 61, "cfc12", cfc12_path, alt = False) #plot the raw values and the rolling mean
plot_figure2(cfc12_data, 62, "cfc12", cfc12_path, alt = False) #plot the seasonal cycle at each station

cfc12_coefficients = find_coefficients(cfc12_data, stations = ["brw"], alt = False) #find coefficients to apply to make the data

plot_figure3(cfc12_data, 63, "cfc12", cfc12_path, alt = False) #plot the residuals (rolling and detrended)
plot_figure4(cfc12_data, cfc12_coefficients, 64, "cfc12", cfc12_path, alt = False) #check that the best fit line actually has 0 intercept
#and is fitting the data well
cfc12_data = do_correction(cfc12_data, cfc12_coefficients, stations = ["brw"]) #apply the offset and ratio
plot_figure5(cfc12_data, cfc12_coefficients, 65, "cfc12", cfc12_path, alt = False) #check that the correction worked
#no calibration ratio applied
plot_figure6(cfc12_data, cfc12_NEEM, 66, "cfc12", cfc12_path, alt = False) #check the overlap with NEEM

#average the dataset together for the final data
cfc12_data, cfc12_egrip = make_final_dataset(cfc12_data, cfc12_NEEM, final_data_directory, "CFC12", alt = False)
#plot the final dataset
plot_figure7(cfc12_egrip, cfc12_NEEM, 67, "cfc12", cfc12_path)

#%% CFC113
cfc113_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\cfc113\\*.txt"
cfc113_files = glob.glob(cfc113_path)
cfc113_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\Final\\cfc113\\"

cfc113_brw_insitu = pd.DataFrame(np.loadtxt(cfc113_files[0], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])

cfc113_sum_insitu = pd.DataFrame(np.loadtxt(cfc113_files[2], usecols = hal_columns), columns = 
                                ["year", "month", "value", "unc", "std"])
cfc113_sum_insitu = clean_data(cfc113_sum_insitu)
cfc113_brw_insitu = clean_data(cfc113_brw_insitu)

cfc113_NEEM = pd.DataFrame(np.loadtxt(cfc113_files[3], usecols = (0,1,2)), columns = 
                                ["date", "value", "std"])

cfc113_discrete = pd.read_table(cfc113_files[1], header = 1)

cfc113_alt_flask, cfc113_brw_flask, cfc113_sum_flask = resample_monthly(cfc113_discrete, "CFC113" )

cfc113_alt_flask = clean_data(cfc113_alt_flask)
cfc113_brw_flask = clean_data(cfc113_brw_flask)
cfc113_sum_flask = clean_data(cfc113_sum_flask)

#make sure that all the new data meshes well with the in situ
plot_figure0(cfc113_brw_flask, cfc113_brw_insitu, 70, "cfc113", "BARROW")
plot_figure0(cfc113_sum_flask, cfc113_sum_insitu, 78, "cfc113", "SUMMIT")

#these don't look like eachother at all and we can see that there is a clear offset, insitu data is ~1ppb larger
#final data is significantly below NEEM so I'll add the mean difference to the flask data

#calculate the offset between the 2 types of data
cfc113_brw_insitu_offset = (cfc113_brw_flask["value"].rolling(12, min_periods = 12, center = True).mean() - 
                        cfc113_brw_insitu["value"].rolling(12, min_periods = 12, center = True).mean()).mean()

cfc113_sum_insitu_offset = (cfc113_sum_flask["value"].rolling(12, min_periods = 12, center = True).mean() - 
                        cfc113_sum_insitu["value"].rolling(12, min_periods = 12, center = True).mean()).mean()

#apply new offset and plot again
cfc113_brw_flask.loc[:,"value"] = cfc113_brw_flask["value"] - cfc113_brw_insitu_offset
cfc113_sum_flask.loc[:, "value"] = cfc113_sum_flask["value"] - cfc113_sum_insitu_offset

#plot to make sure that the data looks good now
plot_figure0(cfc113_brw_flask, cfc113_brw_insitu, 79, "cfc113", "BARROW")
plot_figure0(cfc113_sum_flask, cfc113_sum_insitu, 80, "cfc113", "SUMMIT")


#average each dataset together
cfc113_alt = cfc113_alt_flask.set_index("date")
cfc113_brw = cfc113_brw_flask.set_index("date")
cfc113_sum = cfc113_sum_flask.set_index("date")

cfc113_brw_insitu["flask"] = cfc113_brw["value"]
cfc113_sum_insitu["flask"] = cfc113_sum["value"]

cfc113_brw_insitu["value"] = cfc113_brw_insitu[["value", "flask"]].mean(axis = 1)
cfc113_sum_insitu["value"] = cfc113_sum_insitu[["value", "flask"]].mean(axis = 1)

#make datafile and find all of the values
cfc113_data = make_data_file(cfc113_alt, cfc113_brw_insitu, cfc113_sum_insitu, cfc113_brw_insitu["date"], 
                           cfc113_brw_insitu["month"], cfc113_brw_insitu["year"])
#find the seasonal cycle, subtract it and look at the residuals
cfc113_data = calculate_residuals(cfc113_data, stations = ["alt","brw", "sum"])

plot_figure1(cfc113_data, 71, "cfc113", cfc113_path) #plot the raw values and the rolling mean
plot_figure2(cfc113_data, 72, "cfc113", cfc113_path) #plot the seasonal cycle at each station

cfc113_coefficients = find_coefficients(cfc113_data, stations = ["alt","brw"]) #find coefficients to apply to make the data

plot_figure3(cfc113_data, 73, "cfc113", cfc113_path) #plot the residuals (rolling and detrended)
plot_figure4(cfc113_data, cfc113_coefficients, 74, "cfc113", cfc113_path) #check that the best fit line actually has 0 intercept
#and is fitting the data well
cfc113_data = do_correction(cfc113_data, cfc113_coefficients, stations = ["alt","brw"]) #apply the offset and ratio
plot_figure5(cfc113_data, cfc113_coefficients, 75, "cfc113", cfc113_path) #check that the correction worked
plot_figure6(cfc113_data, cfc113_NEEM, 76, "cfc113", cfc113_path) #check the overlap with NEEM

#average the dataset together for the final data
cfc113_data, cfc113_egrip = make_final_dataset(cfc113_data, cfc113_NEEM, final_data_directory, "CFC113")
#plot the final dataset
plot_figure7(cfc113_egrip, cfc113_NEEM, 77, "cfc113", cfc113_path)

#%% CFC134a
cfc134a_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\cfc134a\\*.txt"
cfc134a_files = glob.glob(cfc134a_path)
cfc134a_path = "C:\\Users\\jrobs\\OneDrive\\Documents\\USRA_2025\\Data\\Final\\cfc134a\\"

raw_cfc134a = pd.read_table(cfc134a_files[0], header = 1)
cfc134a_alt_flask, cfc134a_brw_flask, cfc134a_sum_flask = resample_monthly(raw_cfc134a, "hfc134a")

cfc134a_NEEM = pd.DataFrame(np.loadtxt(cfc134a_files[1]), columns = ["date", "value", "std"])
                    
#clean up the new data
cfc134a_alt_flask = clean_data(cfc134a_alt_flask)
cfc134a_brw_flask = clean_data(cfc134a_brw_flask)
cfc134a_sum_flask = clean_data(cfc134a_sum_flask)

#make datafile and find all of the values
cfc134a_data = make_data_file(cfc134a_alt_flask, cfc134a_brw_flask, cfc134a_sum_flask, 
                              cfc134a_brw_flask["date"], cfc134a_brw_flask["month"], cfc134a_brw_flask["year"])
#find the seasonal cycle, subtract it and look at the residuals
cfc134a_data = calculate_residuals(cfc134a_data)

plot_figure1(cfc134a_data, 81, "cfc134a", cfc134a_path) #plot the raw values and the rolling mean
plot_figure2(cfc134a_data, 82, "cfc134a", cfc134a_path) #plot the seasonal cycle at each station

cfc134a_coefficients = find_coefficients(cfc134a_data) #find coefficients to apply to make the data
plot_figure3(cfc134a_data, 83, "cfc134a", cfc134a_path) #plot the residuals (rolling and detrended)
plot_figure4(cfc134a_data, cfc134a_coefficients, 84, "cfc134a", cfc134a_path) #check that the best fit line actually has 0 intercept
#and is fitting the data well
cfc134a_data = do_correction(cfc134a_data, cfc134a_coefficients) #apply the offset and ratio
plot_figure5(cfc134a_data, cfc134a_coefficients, 85, "cfc134a", cfc134a_path) #check that the correction worked
plot_figure6(cfc134a_data, cfc134a_NEEM, 86, "cfc134a", cfc134a_path) #check the overlap with NEEM

#average the dataset together for the final data
cfc134a_data, cfc134a_egrip = make_final_dataset(cfc134a_data, cfc134a_NEEM, final_data_directory, "CFC134A")
#plot the final dataset
plot_figure7(cfc134a_egrip, cfc134a_NEEM, 87, "cfc134a", cfc134a_path)
          


  