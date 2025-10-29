# -*- coding: utf-8 -*-
"""
Created on Sun Oct 26 16:21:47 2025

@author: jrobs

=== COMPARE MODEL OUTPUTS AT NEEM SITES ===

Comparison of Diffusion Models at different steps of the tuning process. This script compares
the base case model for NEEM with the same model, but with a different calculation for density.
"""
#import libraries
import numpy as np
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

cmap = sns.color_palette(palette='icefire')

#read in datafiles for 1st model
model1_path = "C:\\Users\\jrobs\\eosc448\\FirnAirModel2_fromNEEM\\working_code\\Output_test05_al2_0_uscale01_ust0_aa45\\"
model2_path = "C:\\Users\\jrobs\\eosc448\\FirnAirModel2_fromNEEM\\forwardmodel_8gasfirnairmodel_forNEEM_02\\outputfolder_densityinterp\\"

#%%
#====== 34. Diffusivity Profiles ======
diff1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g0_1000.34"), 
                         columns = ["z", "diff0", "diff"])
diff2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g0_1000.34"), 
                         columns = ["z", "diff0", "diff"])

#make a plot of z against diff and diff1 against diff2
plt.figure(1,clear = True)
plt.plot(diff1["z"], diff1["diff"], label = "NEEM model 1", color = cmap[0])
plt.plot(diff2["z"], diff2["diff"], label = "NEEM model 2", color = cmap[-1], ls= "dashed")
plt.xlabel("z"), plt.ylabel("diffusivity")
plt.legend()

#%%
# ====== 36. Porosity ======
poro_cols =["z","dporo","baro","diff", "ddiff",
        "ddiff*D", "du","-diff(i)*D(n)*grav(n,1)",
        "u*baro", "w","ddiff*D(n)*grav(n,1)","du*baro"]
poro1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g0_1000.36"), 
                         columns = poro_cols)
poro2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g0_1000.36"), 
                         columns = poro_cols)

plt.figure(2, clear = True)
for i, col in enumerate(poro_cols[1:]):
    plt.subplot(3,4,i+1)
    plt.plot(poro1["z"], poro1[col], color = cmap[0], label = "NEEM Model 1")
    plt.plot(poro2["z"], poro2[col], color = cmap[-1], label = "NEEM Model 2", ls = "dashed")
    
    plt.xlabel("z"), plt.ylabel(col)
    plt.title(str(i)+")  "+col)
plt.tight_layout()
plt.legend(bbox_to_anchor = (1.5,0.3))

#%%
# ===== 37. Open Porosity =====
plt.figure(3, clear = True)
openpor_cols = ["z","openpor","dporo","ro", "dro",
         "diff","ddiff", "s", "sc", "w", "u", "du"]
openpor1= pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g0_1000.37"), 
                         columns = openpor_cols)
openpor2= pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g0_1000.37"), 
                         columns = openpor_cols)
for i, col in enumerate(openpor_cols[1:]):
    plt.subplot(3,4,i+1)
    plt.plot(openpor1["z"], openpor1[col], color = cmap[0], label = "NEEM Model 1")
    plt.plot(openpor2["z"], openpor2[col], color = cmap[-1], label = "NEEM Model 2", ls = "dashed")
    
    plt.xlabel("z"), plt.ylabel(col)
    plt.title(str(i)+")  "+col)
plt.tight_layout()
plt.legend(bbox_to_anchor = (1.7,0.3))

#%% 
# ===== 32. Temperature Timeseries for each Gas ====
gases = np.arange(1,9)
temp_columns = ["time", "average_temp_seasonal", "surface_temp", "temp_20m", 
                "temp_50m", "temp_100m", "temp_150m", "average_temp"]
plt.figure(5, clear = True)
cmap = sns.color_palette(palette='tab20')
for gas in gases:
    temp_ts1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.32"), 
                             columns = temp_columns)
    temp_ts2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.32"), 
                             columns = temp_columns)
    plt.subplot(4,2,gas)
    plt.title("Gas #" + str(gas))
    i = 0
    for col in temp_columns[1:]:
        if col == "surface_temp":
            continue
        plt.plot(temp_ts1["time"], temp_ts1[col], label = col, c = cmap[i])
        plt.plot(temp_ts2["time"], temp_ts2[col], c = cmap[i+1], ls = "dashed")
        i +=2 #I'm trying to keep the pairs together
    plt.xlabel("time")
    plt.ylabel("Temperature")
    if gas == 8:
        plt.legend(bbox_to_anchor = (1,1.5))
    plt.tight_layout()
    #make a plot of surface temperatures calculated for each gas
plt.figure(6, clear = True)
for gas in gases:
   temp_ts1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.32"), 
                                 columns = temp_columns)
   temp_ts2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.32"), 
                                 columns = temp_columns)
   plt.subplot(4,2,gas)
   plt.title("Gas #" + str(gas))
   plt.plot(temp_ts1["time"], temp_ts1["surface_temp"], label = col, c = cmap[i])
   plt.plot(temp_ts2["time"], temp_ts2["surface_temp"], c = cmap[i+1], ls = "dashed")
    

#%% 
# ===== 35. Annual Average Temperature Timeseries for each Gas ====
gases = np.arange(1,9)
temp_columns = ["time", "average_temp_seasonal", "average_temp", 
                "temp_50m", "temp_100m", "temp_150m", "surface_temp"]
plt.figure(7, clear = True)
cmap = sns.color_palette(palette='tab20')
for gas in gases:
    temp_ts1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.35"), 
                             columns = temp_columns)
    temp_ts2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.35"), 
                             columns = temp_columns)
    plt.subplot(4,2,gas)
    plt.title("Gas #" + str(gas))
    i = 0
    for col in temp_columns[1:]:
#        if col == "surface_temp":
#            continue
        plt.plot(temp_ts1["time"], temp_ts1[col], label = col, c = cmap[i])
        plt.plot(temp_ts2["time"], temp_ts2[col], c = cmap[i+1], ls = "dashed")
        i +=2 #I'm trying to keep the pairs together
    plt.xlabel("time")
    plt.ylabel("Temperature")
    if gas == 8:
        plt.legend(bbox_to_anchor = (1,1.5))
    plt.tight_layout()
    #make a plot of surface temperatures calculated for each gas
# plt.figure(6, clear = True)
# for gas in gases:
#    temp_ts1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.32"), 
#                                  columns = temp_columns)
#    temp_ts2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.32"), 
#                                  columns = temp_columns)
#    plt.subplot(4,2,gas)
#    plt.title("Gas #" + str(gas))
#    plt.plot(temp_ts1["time"], temp_ts1["surface_temp"], label = col, c = cmap[i])
#    plt.plot(temp_ts2["time"], temp_ts2["surface_temp"], c = cmap[i+1], ls = "dashed")

#%%
# ==== 33. Gas Concentration Profile ====
gasconc_cols = ["z", "gas1", "gas2", "gas3", "gas4", "gas5", "gas6", "gas7", "gas8"]
gasconc1 =  pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g8_1000.33"), 
                         columns = gasconc_cols)
rho1 = openpor1["ro"]
gasconc2= pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g8_1000.33"), 
                         columns = gasconc_cols)
rho2 = openpor2["ro"]

plt.figure(8, clear = True)
c = 0
i = 1
for gas in gasconc_cols[1:]:
    plt.subplot(2,4,i)
    plt.plot(rho1, gasconc1[gas], color = cmap[c])
    plt.plot(rho2, gasconc2[gas], color = cmap[c+1], ls = "dashed")
    plt.xlabel("density")
    plt.ylabel("Gas Concentration")
    plt.title(gas)
    
    i+=1
    c+=2
    
#%% 
# ==== 50+n. Gas concentration Timeseries ====
gases = np.arange(1,9)
gasts_cols = ["time", "surface_conc", "5m_conc", "20m_conc", "lid_conc", "cod_conc", 
              "liz_conc1", "liz_conc2"]
plt.figure(9,clear = True)
for gas in gases:
    gas_ts1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.5"+str(gas)), 
                             columns = gasts_cols)
    gas_ts2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.5"+str(gas)), 
                             columns = gasts_cols)
    
    plt.subplot(2,4,gas)
    i = 0
    for col in gasts_cols[1:]:
        if col == "surface_conc":
            continue
        if col == "5m_conc":
            continue
        if col == "20m_conc":
            continue
        plt.plot(gas_ts1["time"], gas_ts1[col], color = cmap[i], label = col)
        plt.plot(gas_ts2["time"], gas_ts2[col], color = cmap[i+1], ls = "dashed")
        i+=2
    plt.xlabel("time")
    plt.ylabel("Gas Concentration")
    plt.title(gas)
    if gas == 8:
        plt.legend()

#%%
#==== 38. Derivatives =====
deriv_cols = ["z", "d0", "d1", "d2"]
gases = np.arange(1,9)
plt.figure(10, clear = True)
for gas in gases[1:]:
    plt.subplot(2,4,gas)
    derivs1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.38"), 
                             columns = deriv_cols)
    derivs2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.38"), 
                             columns = deriv_cols)
    i = 0
    for col in deriv_cols[1:]:
        plt.plot(derivs1["z"], derivs1[col], color = cmap[i], label = col)
        plt.plot(derivs2["z"], derivs2[col], color = cmap[i+1], ls = "dashed")
        i+=2
    plt.xlabel("density")
    plt.ylabel("Derivative")
    plt.title(gas)
    if gas == 2:
        plt.legend(bbox_to_anchor = (-0.9,1))

#%%
# ==== 39. Depth Temperature and Density ===
gases = np.arange(1,9)
zrot_cols = ["z", "rho", "temp"]
plt.figure(11, clear = True)
ax1 = plt.subplot(1,2,1)
ax2 = plt.subplot(1,2,2)
i = 0
for gas in gases[1:]:
    zrot1 = pd.DataFrame(np.loadtxt(model1_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.39"), 
                             columns = zrot_cols)
    zrot2 = pd.DataFrame(np.loadtxt(model2_path+"NEEMfirnoutput16_test1_g"+str(gas)+"_1000.39"), 
                             columns = zrot_cols)   
    ax1.set_title("temperature profile")
    ax1.plot(zrot1["z"], zrot1["temp"], color = cmap[i], label = gas)
    ax1.plot(zrot2["z"], zrot2["temp"], color = cmap[i+1], ls = "dashed")
    ax1.set_xlabel("depth")
    ax1.set_ylabel("temperature")
    
    ax2.set_title("density profile")
    ax2.plot(zrot1["z"], zrot1["rho"], color = cmap[i], label = gas)
    ax2.plot(zrot2["z"], zrot2["rho"], color = cmap[i+1], ls = "dashed")
    ax2.set_xlabel("depth")
    ax2.set_ylabel("density")  
    i += 2
ax1.legend(bbox_to_anchor = (1,1))