# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:26:14 2024

@author: mcustado
"""

from SALib.analyze.sobol import analyze
from SALib.sample.sobol import sample
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

########################## 3. Run sensitivity analysis ##########################

## Assign input data to variables depending on selected isotope

if iso == 'd18O':
    iso_in = iso_O
    unc_in = unc_O
else:
    iso_in = iso_D
    unc_in = unc_D

## Define the model inputs: hum, temp, dX_P (evap-flux weighted). dX_S, dX_I

problem = {
    'num_vars': 5,
    'names': ['hum', 'temp', 'dX_P', 'dX_S', 'dX_I'],
    'bounds': [[hum-hum_unc, hum+hum_unc],
               [temp-temp_unc, temp+temp_unc],
               [iso_in[3]-unc_in[2], iso_in[3]+unc_in[2]],
               [iso_in[1]-unc_in[3], iso_in[1]+unc_in[3]],
               [iso_in[2]-unc_in[4], iso_in[2]+unc_in[4]]]
    }

## Generate samples

param_values = sample(problem, 1024)

## Run model 

# Initialize empty matrix

Y = np.zeros([param_values.shape[0]])

# Run Sobol method

for i, X in enumerate(param_values):
    print(i)
    alfa = lbf.fractionation_factor_d18O(X[1])  
    ep_eq = (alfa - 1)*1000
    ep_k = lbf.kinetic_en_d18O(X[0])
        
    dX_A = lbf.isotope_atm(X[2], ep_eq, iso_in[0])

    Y[i] = lbf.E_I(X[0], ep_k, ep_eq, alfa, dX_A, X[3], X[4])
    
# Perform analysis

Si = analyze(problem, Y, print_to_console=True)

# ST: Total sensitivity, ST_conf: Confidence interval
# S1: First order sensitivity, S1_conf: Confidence interval
# S2: Second order sensitivity, S2_conf: Confidence interval





