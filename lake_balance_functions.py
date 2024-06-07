# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 21:04:54 2024

@author: mcustado
"""
import lake_balance_functions as lbf

import numpy as np

########## Main function page for lake mass balance analysis #################

# Let:
# h = Humidity
# f = Flux data (discharge/volume change):
# V = Volume of lake

# dX = Isotope value of:
    # L = Lake
    # I = Inflow
    # O = Outlet
    # E = Evaporation
    # P = Precipitation
    # A = Atmosphere (turbulent atmospheric region, based on Craig-Gordon model)
    # S = Isotopic composition of the lake at steady-state
    
# Isotope variable ending in:
    # _O = corresponds to d18O
    # _D = corresponds to dD
    
# ep_k = kinetic enrichment factor
# ep_eq = equilibrium enrichment factor
# alfa = isotopic fractionation factor

# Assumptions: 
    # 1) Lake volume does not change significantly over time
    # 2) Lake is well-mixed
    # 3) Atmospheric conditions constant throughout lake surface

#%% Set up mass balance equations (Gonfiantini, 1981)
    
### Estimate dX_E (isotope value of evaporated water)
# dx_L here can also be dX_O

def isotope_evap (h, ep_k, ep_eq, alfa, dX_A, dX_L):
    dX_E = (((dX_L-ep_eq)/alfa)-(h*dX_A)-ep_k)/(1-h+(0.001*ep_k))
    return dX_E

# Estimate isotopic composition of atm moisture
# Atmospheric moisture upwind of the lake is assumed to be in equilibrium with mean annual precipitation or if site is seasonal, precipitation during the evaporation season. Ideally, evaporation flux-weighted precipitation isotope data (see also Gibson et al., 2008, 2015)
# k parameter = seasonality factor. 0.5 (highly seasonal) to 1 (non-seasonal). Gibston, et al., 2015

def isotope_atm (precip, ep_eq, k): 
    dX_A = (precip-(k*ep_eq))/(1+(0.001*k*ep_eq))
    return dX_A

### Estimate isotopic composition of lake at steady state (dX_S)
# dX_L here can also be dX_O

def mass_balance_ss2 (h, ep_k, ep_eq, alfa, dX_A, dX_S, dX_I):
    a = ((h*dX_A)+ep_k+(ep_eq/alfa))/(1-h+(0.001*ep_k))
    b = (h-0.001*(ep_k+(ep_eq/alfa)))/(1-h+(0.001*ep_k))
    x = (dX_S-dX_I)/(a-(b*dX_S))
    dX_LS = ((x*a)+dX_I)/(1+(b*x))
    return dX_LS, x, a, b # Returns lake steady state isotopic composition, X, terms A and B (See Equations 5 and 6 in Custado, et al. 2024)

def E_I (h, ep_k, ep_eq, alfa, dX_A, dX_S, dX_I):
    a = ((h*dX_A)+ep_k+(ep_eq/alfa))/(1-h+(0.001*ep_k))
    b = (h-0.001*(ep_k+(ep_eq/alfa)))/(1-h+(0.001*ep_k))
    E_I = (dX_S-dX_I)/(a-(b*dX_S))
    return E_I # Returns just X

def mass_balance_ssx (h, ep_k, ep_eq, alfa, dX_A, dX_I, x):
    a = ((h*dX_A)+ep_k+(ep_eq/alfa))/(1-h+(0.001*ep_k))
    b = (h-0.001*(ep_k+(ep_eq/alfa)))/(1-h+(0.001*ep_k))
    dX_LS = ((x*a)+dX_I)/(1+(b*x))
    limit = a/b
    return dX_LS, limit # Returns lake steady state isotopic composition and theoretical maximum enrichment of lake

### Construct hydrological balance by simultaneously calculating for evaporation, groundwater discharge, X, humidity, and isotpic composition of inflow

def hydro_balance(vars, dX_S_O, dX_S_D, dX_A_O, dX_A_D, f_inlet, f_creek, f_precip, f_outlet, alfa_O, ep_eq_O, alfa_D, ep_eq_D):
    
    # define constants
    dX_inlet_O = -16.5680552351257
    dX_creek_O = -16.6427693908244
    dX_precip_O = -14.5993690452293
    dX_gwater_O = -17.7983116883116
    dX_inlet_D = -125.869415352418
    dX_creek_D = -126.203603690239
    dX_precip_D = -105.704124214273
    dX_gwater_D = -135.735649350649

    # define unknowns
    
    f_gwater, f_evap, x_, h, dX_I_O, dX_I_D = vars
    
    eq1 = f_inlet + f_creek + f_precip + f_gwater - f_outlet - f_evap
    eq2 = x_*(f_inlet + f_creek + f_gwater + f_precip) - f_evap
    
    eq3 = ((dX_inlet_O*f_inlet + dX_creek_O*f_creek + dX_precip_O*f_precip + dX_gwater_O*f_gwater) / (f_inlet + f_creek + f_precip + f_gwater)) - dX_I_O
    eq4 = ((dX_inlet_D*f_inlet + dX_creek_D*f_creek + dX_precip_D*f_precip + dX_gwater_D*f_gwater) / (f_inlet + f_creek + f_precip + f_gwater)) - dX_I_D
    
    ep_k_O = 14.2*(1-h)
    ep_k_D = 12.5*(1-h)

    eq5 = (((dX_S_O - dX_I_O)*(1-h+(0.001*ep_k_O))) / (h*(dX_A_O - dX_S_O) + (ep_k_O + (ep_eq_O/alfa_O))*(0.001*dX_S_O + 1))) - x_
    eq6 = (((dX_S_D - dX_I_D)*(1-h+(0.001*ep_k_D))) / (h*(dX_A_D - dX_S_D) + (ep_k_D + (ep_eq_D/alfa_D))*(0.001*dX_S_D + 1))) - x_
    
    return np.array([eq1, eq2, eq3, eq4, eq5, eq6])

### Function for simulation of X used in Section 6.3 (Custado, et al. 2024)

def calc_x(vars, dX_S_O, dX_I_O, dX_P_O, h, temp_): #dX_E_O, dX_E_D, 
    x_ = vars
    
    alfa_O = lbf.fractionation_factor_d18O(temp_)
    ep_eq_O = (alfa_O - 1)*1000

    ep_k_O = 14.2*(1-h)

    dX_A_O = lbf.isotope_atm(dX_P_O, ep_eq_O, 1)

    eq1 = (((dX_S_O - dX_I_O)*(1-h+(0.001*ep_k_O))) / (h*(dX_A_O - dX_S_O) + (ep_k_O + (ep_eq_O/alfa_O))*(0.001*dX_S_O + 1))) - x_

    return np.array(eq1)


#%% Set up equations for fractionation and enrichment factors

def fractionation_factor_d18O (temp): # Temperature input in deg_C
    alpha = np.exp((-7.685/(10**3)) + (6.7123/(273.15 + temp)) - (1666.4/((273.15 + temp)**2)) + (350410/((273.15 + temp)**3)))
    return alpha

def fractionation_factor_dD (temp): # Temperature input in deg_C
    alpha = np.exp((1158.8*(((273.15 + temp)**3)/(10**12))) - (1620.1*(((273.15 + temp)**2)/(10**9))) + (794.84*((273.15 + temp)/(10**6))) - (161.04/(10**3)) + (2999200/((273.15 + temp)**3)))
    return alpha

def kinetic_en_d18O(humidity):
    ep_k = 14.2*(1-humidity)
    return ep_k

def kinetic_en_dD(humidity):
    ep_k = 12.5*(1-humidity)
    return ep_k

#%% Print data functions

# Data for first mass balance calculations

def print_results_calcs1(iso, hum, temp, ep, ep_k, alpha, ds, dX_A, dX_E, xs):
    
    print ("Inputs to model:",
           "\nIsotope:\t",iso,
           "\nHumidity:\t",hum,
           "\nTemperature (deg C):\t",temp,
           "\nEq. enrichment factor:\t",ep,
           "\nKinetic enrichment factor:\t",ep_k,
           "\nEq. fractionation factor:\t",alpha,
           "\nIsotope - atmosphere:\t",dX_A,
           "\nIsotope - lake, steady-state (available lake data):\t",ds[2],
           "\nIsotope - inflow:\t",ds[0],
           "\n\nOutput:",
           "\nIsotope - evaporation (calculated):\t",dX_E,
           "\nX:\t",xs)
    
# Data for uncertainty calculations

def print_results_unc(iso, sim, x_dist, dX_E_dist, dX_A_dist, E_dist):
        
    print("Isotope:\t",iso,
        "\nNumber of simulations:\t", sim)
    
    print("\nX output:",
        "\nmean output:\t", np.nanmean(x_dist),
        "\nmedian output:\t", np.median(x_dist),
        "\nstdev output:\t", np.nanstd(x_dist),
        "\n15.9 perc output:\t", np.percentile(x_dist, 15.9),
        "\n84.1 perc output:\t", np.percentile(x_dist, 84.1),
        "\nminimum output:\t", np.min(x_dist),
        "\nmaximum output:\t", np.max(x_dist))
    
    print("\ndX_E output:",
        "\nmean output:\t", np.nanmean(dX_E_dist),
        "\nmedian output:\t", np.median(dX_E_dist),
        "\nstd output:\t", np.nanstd(dX_E_dist),
        "\n15.9 perc output:\t", np.percentile(dX_E_dist, 15.9),
        "\n84.1 perc output:\t", np.percentile(dX_E_dist, 84.1),
        "\nminimum output:\t", np.min(dX_E_dist),
        "\nmaximum output:\t", np.max(dX_E_dist))
        
    print("\ndX_A output:",
        "\nmean output:\t", np.nanmean(dX_A_dist),
        "\nmedian output:\t", np.median(dX_A_dist),
        "\nstd output:\t", np.nanstd(dX_A_dist),
        "\n15.9 perc output:\t", np.percentile(dX_A_dist, 15.9),
        "\n84.1 perc output:\t", np.percentile(dX_A_dist, 84.1),
        "\nminimum output:\t", np.min(dX_A_dist),
        "\nmaximum output:\t", np.max(dX_A_dist))
    
    print("\nE output:",
        "\nmean output:\t", np.nanmean(E_dist),
        "\nmedian output:\t", np.median(E_dist),
        "\nstd output:\t", np.nanstd(E_dist),
        "\n15.9 perc output:\t", np.percentile(E_dist, 15.9),
        "\n84.1 perc output:\t", np.percentile(E_dist, 84.1),
        "\nminimum output:\t", np.min(E_dist),
        "\nmaximum output:\t", np.max(E_dist))
    