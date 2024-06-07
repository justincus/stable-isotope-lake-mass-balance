# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 19:03:05 2024

@author: mcustado
"""

import scipy.optimize as opt
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

atm_O = lbf.isotope_atm(precip_O, ep_O, 1) # k=1 seasonality constant
atm_D = lbf.isotope_atm(precip_D, ep_D, 1) # k=1 seasonality constant

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

### Back calculate f_gwater, f_evap, x_, h, dX_I_O, dX_I_D  ###

# Input known discharge values (m3/yr):
    
f_inlet2 = 317867056.26687
f_creek2 = 145049044.909472
f_precip2 = 107856450.331302
f_outlet2 = 352639058.615613

# Input initial guesses for output variables (f_gwater, f_evap, x_, h, dX_I_O, dX_I_D)

initial_guesses = [0, 231211349.583575, 0.435, 0.76, -16.2150279446561, -122.143792918018]

# Run solver

roots = opt.fsolve(lbf.hydro_balance, initial_guesses, args = (lake_O, lake_D, atm_O, atm_D, f_inlet2, f_creek2, f_precip2, f_outlet2, alpha_O, ep_O, alpha_D, ep_D)) #, method='hybr')

f_gwater_soln, f_evap_soln, x_soln, h_soln, dX_I_O_soln, dX_I_D_soln = roots[0], roots[1], roots[2], roots[3], roots[4], roots[5]

print("f_gwater_soln:", f_gwater_soln, "\nf_evap_soln:", f_evap_soln, 
      "\nx_soln:", x_soln, "\nh_soln:", h_soln,
      "\ndX_I_O_soln:", dX_I_O_soln, "\ndX_I_D_soln:", dX_I_D_soln,
      
      "\n\nf_inlet:", f_inlet2, "\nf_creek:", f_creek2,
      "\nf_precip:", f_precip2, "\nf_outlet:", f_outlet2)
