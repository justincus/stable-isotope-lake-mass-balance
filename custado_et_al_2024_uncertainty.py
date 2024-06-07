# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:26:14 2024

@author: mcustado
"""

import numpy as np
import lake_balance_functions as lbf

########################## 1. Input parameters ##########################

## Choose which isotope to analyze

iso = 'dD' # Select stable isotope for analysis (dD or d18O)

## Climate data

hum = 0.62
temp = 11.15

## Isotope data

# d18O 

k_O = 1

dX_P_O = -11.70
dX_S_O = -8.75978345841666
dX_I_O = -16.2152393388515

iso_O = [k_O, dX_S_O, dX_I_O, dX_P_O] # Place in one array

# dD 

k_D = 1

dX_P_D = -84.02
dX_S_D = -86.4422222222222
dX_I_D = -122.145607652468

iso_D = [k_D, dX_S_D, dX_I_D, dX_P_D] # Place in one array


########################## 2. Input uncertainties ##########################

## Input uncertanties

# Climate data

hum_unc = hum*0.05 #mean*(percent/100) # Uncertainty for humidity (%)
temp_unc = 0.2 # Uncertainty for temperature (degC)

# d18O

dX_P_O_unc = abs(dX_P_O*0.003) # Uncertainty for isotopic composition of evaporation flux-weighted precipitation (per mil)
dX_S_O_unc = 0.1  # Uncertainty for isotopic composition of steady state lake (per mil)
dX_I_O_unc = 0.0454711273463026 # Uncertainty for isotopic composition of total inflow (per mil)

unc_O = [hum_unc, temp_unc, dX_P_O_unc, dX_S_O_unc, dX_I_O_unc] # Place in one array

# dD

dX_P_D_unc = abs(dX_P_D*0.01) # Uncertainty for isotopic composition of evaporation flux-weighted precipitation (per mil)
dX_S_D_unc = 0.5 # Uncertainty for isotopic composition of steady state lake (per mil)
dX_I_D_unc = 0.338982073484539 # Uncertainty for isotopic composition of total inflow (per mil)

unc_D = [hum_unc, temp_unc, dX_P_D_unc, dX_S_D_unc, dX_I_D_unc] # Place in one array

########################## 3. Calculate combined uncertainty ##########################

## Assign input data to variables depending on selected isotope

if iso == 'd18O':
    iso_in = iso_O
    unc_in = unc_O
else:
    iso_in = iso_D
    unc_in = unc_D

## Input parameters: hum, temp, dX_P (evap-flux weighted). dX_S, dX_I

# Input # of simulations

sim = 100000

## Initialize arrays of input distributions

hum_dist = np.random.default_rng().uniform(hum-unc_in[0], hum+unc_in[0], sim)
temp_dist = np.random.default_rng().uniform(temp-unc_in[1], temp+unc_in[1], sim)

dX_P_dist = np.random.default_rng().uniform(iso_in[3]-unc_in[2], iso_in[3]+unc_in[2], sim)
dX_S_dist = np.random.default_rng().uniform(iso_in[1]-unc_in[3], iso_in[1]+unc_in[3], sim)
dX_I_dist = np.random.default_rng().uniform(iso_in[2]-unc_in[4], iso_in[2]+unc_in[4], sim)

## Initialize output arrays

x_dist = []
dX_A_dist = []
dX_E_dist = []
E_dist = []

## Run mass balance function simulations

for ind in np.arange(0,sim,1):
    print(ind)

    # Calculate fractionation factors

    alfa = lbf.fractionation_factor_d18O(temp_dist[ind])  
    ep_eq = (alfa - 1)*1000
    ep_k = lbf.kinetic_en_d18O(hum_dist[ind])
    
    # Calculate isotopic composition of atmospheric moisture
        
    dX_A = lbf.isotope_atm(dX_P_dist[ind], ep_eq, iso_in[0])
    dX_A_dist.append(dX_A)
    
    # Calculate isotopic composition of evaporate
    
    dX_E = lbf.isotope_evap(hum_dist[ind], ep_k, ep_eq, alfa, dX_A, dX_S_dist[ind])
    dX_E_dist.append(dX_E)

    # Calculate X

    x_ = lbf.E_I(hum_dist[ind], ep_k, ep_eq, alfa, dX_A, dX_S_dist[ind], dX_I_dist[ind]) 

    x_dist.append(x_)
    
    # Calculate actual evaporation rate (m3/day)
    
    E_ = x_*570772551.507645 # Input total annual volumetric inflow (m3/day)
    
    E_dist.append(E_)
    
## Print results

lbf.print_results_unc(iso, sim, x_dist, dX_E_dist, dX_A_dist, E_dist)



