# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:00:30 2024

@author: mcustado
"""
import lake_balance_functions as lbf

################ 1. Input parameters ################

# Choose which isotope to analyze

iso = 'dD' # Select stable isotope for analysis (dD or d18O)

# Climate data
temp = 11.15 # Input evaporation-flux weighted temperature
hum = 0.62 # Input evaporation-flux weighted humidity

# d18O data
influx_O = -16.2152393388515
lake_O = -8.75978345841666
precip_O = -11.70 # Evaporation flux-weighted

# dD isotope data
influx_D = -122.145607652468
lake_D = -86.4422222222222
precip_D = -84.02 # Evaporation flux-weighted

# Calculate equilibrium fractionation (alpha) and enrichment (ep) factors

alpha_O = lbf.fractionation_factor_d18O(temp) 
alpha_D = lbf.fractionation_factor_dD(temp) 

ep_O = (alpha_O - 1)*1000
ep_D = (alpha_D - 1)*1000

# Calculate kinetic enrichment factor (ep_k)

ep_k_O = lbf.kinetic_en_d18O(hum)
ep_k_D = lbf.kinetic_en_dD(hum)

# Calculate dX_A and dX_E

atm_O = lbf.isotope_atm(precip_O, ep_O, 1) # k seasonality constant
atm_D = lbf.isotope_atm(precip_D, ep_D, 1) # k seasonality constant

evap_O = lbf.isotope_evap(hum, ep_k_O, ep_O, alpha_O, atm_O, lake_O) 
evap_D = lbf.isotope_evap(hum, ep_k_D, ep_D, alpha_D, atm_D, lake_D) 

# Assign input data to variables depending on selected isotope

if iso == 'd18O':
    ds = [influx_O, lake_O, precip_O, ep_k_O, ep_O, alpha_O, atm_O, evap_O]
else:
    ds = [influx_D, lake_D, precip_D, ep_k_D, ep_D, alpha_D, atm_D, evap_D]

influx = ds[0]
lake = ds[1]
precip = ds[2]
ep_k = ds[3]
ep = ds[4]
alpha = ds[5]
atm = ds[6]
evap = ds[7]

################ 2. Run mass balance calculations ################

# Run mass balance equation

dX_LS, xs, a, b = lbf.mass_balance_ss2(hum, ep_k, ep, alpha, atm, lake, influx)

# Print results

lbf.print_results_calcs1(iso, hum, temp, ep, ep_k, alpha, ds, atm, evap, xs)

