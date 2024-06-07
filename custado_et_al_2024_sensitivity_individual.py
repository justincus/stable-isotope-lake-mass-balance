# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:26:14 2024

@author: mcustado
"""

import numpy as np
import matplotlib.pyplot as plt
import lake_balance_functions as lbf

########################## 1. Input parameters ##########################

## Choose which isotope to analyze

iso = 'dD' # Select stable isotope for analysis (dD or d18O)

# Provide climate data

hum = 0.62
temp = 11.15
precip_O = -11.70 # evap flux-weighted
precip_D = -84.02 # evap flux-weighted

# Provide fractionation factors data

alfa_O = lbf.fractionation_factor_d18O(temp)  
ep_eq_O = (alfa_O - 1)*1000
ep_k_O = lbf.kinetic_en_d18O(hum)

frac_O = [alfa_O, ep_eq_O, ep_k_O] # d18O

alfa_D = lbf.fractionation_factor_dD(temp) 
ep_eq_D = (alfa_D - 1)*1000
ep_k_D = lbf.kinetic_en_dD(hum)

frac_D = [alfa_D, ep_eq_D, ep_k_D] # dD

# Provide input isotopic composition data

k_O = 1
dX_A_O = lbf.isotope_atm(precip_O, ep_eq_O, k_O)
dX_S_O = -8.75978345841666
dX_I_O = -16.2152393388515
dX_P_O = -11.70

iso_O = [k_O, dX_A_O, dX_S_O, dX_I_O, dX_P_O]

k_D = 1
dX_A_D = lbf.isotope_atm(precip_D, ep_eq_D, k_D)
dX_S_D = -86.4422222222222
dX_I_D = -122.145607652468
dX_P_D = -84.02

iso_D = [k_D, dX_A_D, dX_S_D, dX_I_D, dX_P_D]

# Provide input uncertanties

hum_unc = hum*0.05 #mean*percent/100
dX_P_O_unc = abs(dX_P_O*0.003) #mean*percent/100 #1.44029825 # 4.8 #stdev 
dX_P_D_unc = abs(dX_P_D*0.01) #mean*percent/10010.90144133 # 36.9 #stdev 

temp_unc = 0.2 #degC
dX_S_O_unc = 0.1 #0.025 #per mil #0.726266967 #stdev
dX_I_O_unc = 0.0454711273463026 #0.0150159 #combined_precision #1.44029825 #stdev
dX_S_D_unc = 0.5 # 0.1 #per mil #3.422022721 #stdev 0.1 #per mil
dX_I_D_unc = 0.338982073484539 #combined_precision per mil #10.90144133 #stdev 0.1 #per mil

unc_O = [hum_unc, temp_unc, dX_P_O_unc, dX_S_O_unc, dX_I_O_unc]
unc_D = [hum_unc, temp_unc, dX_P_D_unc, dX_S_D_unc, dX_I_D_unc]

########################## 2. Run simulations ##########################

## Assign input data to variables depending on selected isotope

if iso == 'd18O':
    iso_in = iso_O
    frac_in = frac_O
    unc_in = unc_O
else:
    iso_in = iso_D
    frac_in = frac_D
    unc_in = unc_D

## Humidity vs X 
## Vary humidity from 0.5 to 1

# Input # of simulations

sim = 1000000

# Initialize array of input distribution

hum_dist = np.random.default_rng().uniform(0.5, 1, sim)

# Initialize output array

x_dist = []

# Run simulations

for ind in np.arange(0,sim,1):
    print(ind)

    x_ = lbf.E_I(hum_dist[ind], frac_in[2], frac_in[1], frac_in[0], iso_in[1], iso_in[2], iso_in[3]) #[alfa_O, ep_eq_O, ep_k_O] iso_O = [k_O, dX_A_O, dX_S_O, dX_I_O]
    x_dist.append(x_)
    
# Plot

input_x = hum_dist
output_y = x_dist
input_variable = "hum_dist"
output_variable = "X"
fname = iso+"_"+input_variable+"_"+output_variable

fig, ax1 = plt.subplots(figsize=(17, 10))

# Creating a twin axis sharing the x-axis
ax2 = ax1.twinx()

# Plotting on the first axis
ax1.plot(output_y, input_x, 'o', zorder=1)  # Ensure line is plotted on top
ax1.set_ylabel('Input: Humidity', fontsize=25)
ax1.set_xlabel('Output: X (Evaporation/Inflow)', fontsize=25)
ax1.tick_params(axis='both', which='major', labelsize=25)

# Plotting histogram on the twin axis
weights = np.ones_like(output_y) / len(output_y)
ax2.hist(output_y, 50, weights=weights, color='gray', zorder=2, alpha=0.5)
ax2.set_ylabel('Probability', fontsize=25)
ax2.tick_params(axis='both', which='major', labelsize=25)

ax1.set_zorder(1)  # default zorder is 0 for ax1 and ax2
ax1.set_frame_on(False)  # prevents ax1 from hiding ax2

# Title
plt.suptitle(iso + " Input: " + input_variable + " Output: " + output_variable, fontsize=20, y=1.02)
plt.tight_layout()

plt.show()

## dX_A vs X 
## Vary dX_A by varying input evaporation flux-weighted dX_P

# Input # of simulations

sim = 100000

# Initialize array of input distribution

precip_in = np.random.default_rng().uniform(iso_in[4]*1.2, iso_in[4]*0.8, sim)
dX_A_dist = []

# Initialize output array

x_dist = []

# Run simulations

for ind in np.arange(0,sim,1):
    print(ind)

    dX_A_dist_ = lbf.isotope_atm(precip_in[ind], frac_in[1], 1)
    dX_A_dist = np.append(dX_A_dist, dX_A_dist_)

    x_ = lbf.E_I(hum, frac_in[2], frac_in[1], frac_in[0], dX_A_dist_, iso_in[2], iso_in[3]) #[alfa_O, ep_eq_O, ep_k_O] iso_O = [k_O, dX_A_O, dX_S_O, dX_I_O]
    x_dist.append(x_)
    
# Plot

input_x = dX_A_dist
output_y = x_dist
input_variable = "dX_A_dist"
output_variable = "X"
fname = iso+"_"+input_variable+"_"+output_variable

fig, ax1 = plt.subplots(figsize=(17, 10))

# Creating a twin axis sharing the x-axis
ax2 = ax1.twinx()

# Plotting on the first axis
ax1.plot(output_y, input_x, 'o', zorder=1)  # Ensure line is plotted on top
ax1.set_ylabel('Input: δA, ‰', fontsize=25)
ax1.set_xlabel('Output: X (Evaporation/Inflow)', fontsize=25)
ax1.tick_params(axis='both', which='major', labelsize=25)

# Plotting histogram on the twin axis
weights = np.ones_like(output_y) / len(output_y)
ax2.hist(output_y, 50, weights=weights, color='gray', zorder=2, alpha=0.5)
ax2.set_ylabel('Probability', fontsize=25)
ax2.tick_params(axis='both', which='major', labelsize=25)

ax1.set_zorder(1)  # default zorder is 0 for ax1 and ax2
ax1.set_frame_on(False)  # prevents ax1 from hiding ax2

# Title
plt.suptitle(iso + " Input: " + input_variable + " Output: " + output_variable, fontsize=20, y=1.02)
plt.tight_layout()

plt.show()




