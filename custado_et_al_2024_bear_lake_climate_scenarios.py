# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 19:03:05 2024

@author: mcustado
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import lake_balance_functions as lbf

################ 1. Input parameters ################

## Provide climate data

hum = 0.62 #current
temp = 11.15 #current

period = ['lig','current','future','glacial']
# all arrays = [LIG, current, future, glacial (LGP)]

## Provide input isotopic composition data

k = [1,1,1,1] # seasonality constant for lig, current, future, and glacial periods
dX_P = [-11.7, -11.7, -11.7, -11.7] # isotopic composition of precipitation for lig, current, future, and glacial periods
dX_S = [-7.21, -8.76, -8.76, -13.13] # steady-staete isotopic composition of bear lake for lig, current, future, and glacial periods
dX_I = [-15.77,-16.22,-16.22,-16.22] # isotopic composition of total inflow for lig, current, future, and glacial periods
E_I = 0.38 # calculated modern evaporation / inflow (X)

## Provide humidity and temperature changes for each scenario

hum_dec = [0.1, 0, 0.1, 0.1] # humidity decrease for lig, current, future, and glacial periods
temp_inc = [1, 0, 1, -6] # temperature decrease (oC) for lig, current, future, and glacial periods

## Provide input uncertanties

hum_unc = 0.03
temp_unc = 0.3
dX_P_O_unc = 2 # per mil
dX_S_O_unc = [2,1,1,0.9] # per mil
dX_I_O_unc = 2 # per mil
E_I_unc = 0.1

################ 2. Run simulations ################

## Input number of simulations

sim = 100000

## Initialize input distributions

# initialize index
index = np.arange(0,4,1) 

# initialize input arrays
hum_in = []
temp_in = []
lake_in = []

# input input distributions: loop through different periods, create temp array for each input parameter (hum, temp, lake)
for i in index:
    print(i)

    hum_temp = np.random.default_rng().normal(hum-hum_dec[i], hum_unc, sim)
    hum_in.append(hum_temp)

    temp_temp = np.random.default_rng().normal(temp+temp_inc[i], temp_unc, sim)
    temp_in.append(temp_temp)
    
    lake_temp = np.random.default_rng().uniform(dX_S[i]-dX_S_O_unc[i], dX_S[i]+dX_S_O_unc[i], sim)
    lake_in.append(lake_temp)

## Calculate for x

# initialize output array
x_array = []

for i in index:
    x_temp = [] # create temporary output array
    
    for ind in np.arange(0,sim,1):
        print(ind)
        roots = opt.fsolve(lbf.calc_x, E_I, args = (lake_in[i][ind], dX_I[i], dX_P[i], hum_in[i][ind], temp_in[i][ind]))
        x_temp = np.append(x_temp, roots[0])
        
    x_array.append(x_temp)

# print results for X in each period

for i in index:

    print("period: ", period[i])
    print("number of simulations: ", sim)
    print("output: x")
    
    print("mean output: \t", np.nanmean(x_array[i]))
    print("15.9 perc output: \t", np.percentile(x_array[i], 15.9))
    print("84.1 perc output: \t", np.percentile(x_array[i], 84.1))
    print("minimum output: \t", np.min(x_array[i]))
    print("maximum output: \t", np.max(x_array[i]), "\n")

################ 3. Plot dX_S, humidity, temperatuve vs X in different scenarios ################

for i in index: # loop through the four scenarios (first to last output plots: LIG, current, future, and glacial (LGP) scenarios)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(24,7), sharey=True)
    
    ax1.scatter(lake_in[i], x_array[i], alpha=0.1, color='#c0c0c0')
    ax1.set_ylabel('X (Evaporation/Inflow)', fontsize=25)
    ax1.set_xlabel('Lake δ$^1$$^8$O (‰)', fontsize=25)
    m,b = np.polyfit(lake_in[i], x_array[i], deg=1)
    y = m*lake_in[i] + b
    ax1.plot(lake_in[i], y, linewidth=2, zorder=1)
    ax1.scatter(dX_S[i], np.mean(x_array[i]), color='black', s=200)
    ax1.tick_params(axis='x', labelsize=25)
    ax1.tick_params(axis='y', labelsize=25)
    ax1.grid()
    
    ax2.scatter(hum_in[i], x_array[i], alpha=0.1, color='#c0c0c0')
    ax2.set_xlabel('Humidity', fontsize=25)
    ax2.tick_params(axis='x', labelsize=25)
    ax2.grid()
    
    ax3.scatter(temp_in[i], x_array[i], alpha=0.1, color='#c0c0c0')
    ax3.set_xlabel('Temperature (ºC)', fontsize=25)
    ax3.tick_params(axis='x', labelsize=25)
    ax3.grid()

    print(period[i])
    print("slope: \t",m)
    print("y-int: \t",b)
    
plt.tight_layout()
# plt.savefig('.....\\'+period[i]+'.png', bbox_inches="tight", dpi=600)

################ 4. Plot all dX_S, humidity, temperatuve vs X scenarios in one field ################

colors = ['#09A603', '#D9B504', '#D90404', '#0583F2']

legend = ["LIG", "Current", "Future", "LGP"]

fig, ax1 = plt.subplots(figsize=(8,7))

for i in index: # loop through the four scenarios (first to last output plots: LIG, current, future, and glacial (LGP) scenarios)
    
    ax1.scatter(lake_in[i], x_array[i], alpha=0.1, color=colors[i], label = legend[i])
    ax1.set_ylabel('X (Evaporation/Inflow)', fontsize=15)
    ax1.set_xlabel('Lake δ$^1$$^8$O (‰)', fontsize=15)
    m,b = np.polyfit(lake_in[i], x_array[i], deg=1)
    y = m*lake_in[i] + b
    ax1.plot(lake_in[i], y, linewidth=2, zorder=1)
    ax1.scatter(dX_S[i], np.mean(x_array[i]), color='black', s=200)
    ax1.tick_params(axis='x', labelsize=15)
    ax1.tick_params(axis='y', labelsize=15)
    leg = ax1.legend(fontsize=15, markerscale = 2)
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
    ax1.grid(visible=True, alpha = 0.5)
    
plt.tight_layout()
# plt.savefig('.....\\all_scenarios_plot.png', bbox_inches="tight", dpi=600)

fig, (ax2, ax3) = plt.subplots(1,2, figsize=(16,7), sharey=True)

for i in index: # loop through the four scenarios (first to last output plots: LIG, current, future, and glacial (LGP) scenarios)
        
    ax2.scatter(hum_in[i], x_array[i], alpha=0.1, color=colors[i], label = legend[i])
    ax2.set_ylabel('X (Evaporation/Inflow)', fontsize=15)
    ax2.set_xlabel('Humidity', fontsize=20)
    ax2.tick_params(axis='y', labelsize=20)
    leg = ax2.legend(fontsize=15, markerscale = 2, loc = 'upper left')
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
    ax2.tick_params(axis='x', labelsize=20)
    ax2.grid(visible=True, alpha = 0.5)
    
    ax3.scatter(temp_in[i], x_array[i], alpha=0.1, color=colors[i], label = legend[i])
    ax3.set_xlabel('Temperature (ºC)', fontsize=20)
    ax3.tick_params(axis='x', labelsize=20)
    ax3.grid(visible=True, alpha = 0.5)
    
plt.tight_layout()
# plt.savefig('.....\\all_scenarios_plots_for_supp2.png', bbox_inches="tight", dpi=600)
