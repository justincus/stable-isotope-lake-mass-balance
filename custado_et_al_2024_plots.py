# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:06:44 2024

@author: mcustado
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import lake_balance_functions as lbf
import geopandas as gp

################ Load masterlist of data and relevant shapefiles ################

bl = pd.read_csv(r'.......\\BL_master_list.csv') # Master isotope data list [included in datasets provided]
bl['Sample_Collection_Date'] = pd.to_datetime(bl['Sample_Collection_Date'])
df = gp.read_file('..........\\Clipped Bear Lake shapefile\\bl_clipped.shp') # Clipped shapefile of Bear Lake [included in datasets provided]
rv = gp.read_file(r'..........\\\Lakes_and_Rivers_Shapefile_NA_Lakes_and_Rivers_data_hydrography_l_rivers_v2\Lakes_and_Rivers_Shapefile\NA_Lakes_and_Rivers\data\hydrography_l_rivers_v2.shp') ## Source: https://www.sciencebase.gov/catalog/item/4fb55df0e4b04cb937751e02
rv_wgs84 = rv.to_crs({'init': 'epsg:4326'}) 
basin = gp.read_file('..........\\\Great Basin\Shape\WBDHU12.shp') ## Source: https://www.sciencebase.gov/catalog/item/52c7d4cbe4b0a753c7d3c586
stations_used = pd.read_csv(r'..........\\\stations_used.csv') # Hydrological stations used in the calculations [included in datasets provided]

################ Figure 1: Sample map ################

# initialize an axis
plt.rcParams.update({'font.size': 20})
fig, ax = plt.subplots(figsize=(10,10))
# plot map on axis
df.plot(alpha=0.05, edgecolor='k', color='lightgrey', ax=ax)
ax.set_ylim(40.5, 43)
ax.set_xlim(-112.5, -110.5)

rv_wgs84.loc[rv_wgs84['NAMEEN'].str.contains(r'Bear',na=False)].plot(color = '#86BBD8', ax=ax)

basin.loc[basin['name']=="Bear Lake"].plot(color='#86BBD8', edgecolor = 'black', linewidth=2, ax=ax)

bl.loc[(bl['Type']=="Snow_pit")].plot(x="Lon", y="Lat", kind="scatter", s = 150, color = '#05D5FA', label = 'Snow pit', linewidth = 0.5, edgecolor='black', ax=ax)
bl.loc[(bl['Type']=="Ground") | (bl['Type']=="Spring")].plot(x="Lon", y="Lat", kind="scatter", s = 150, color = '#C4A484', label = 'Ground and spring', linewidth = 0.5, edgecolor='black', ax=ax)
bl.loc[(bl['Type']=="Canal")].plot(x="Lon", y="Lat", kind="scatter", s = 150, color = '#9EE493', label = 'Canal', linewidth = 0.5, edgecolor='black', ax=ax)
bl.loc[(bl['Type']=="River_or_stream")].plot(x="Lon", y="Lat", kind="scatter", s = 150, color = '#3D7BBA', label = 'River/stream', linewidth = 0.5, edgecolor='black', ax=ax)
bl.loc[(bl['Type']=="Lake")].plot(x="Lon", y="Lat", marker = 'D', kind="scatter", s = 150, color = '#3D7BBA', label = 'Lake', linewidth = 0.5, edgecolor='black', ax=ax)
bl.loc[(bl['Type']=="Precipitation")].plot(x="Lon", y="Lat",kind="scatter", s = 150, color = '#DA70D6', label = 'Precipitation', linewidth = 0.5, edgecolor='black', ax=ax)
bl.loc[(bl['Data_source']=="Project")].plot(x="Lon", y="Lat",kind="scatter", s = 35, color = '#FF8A00', label = '2022/2023 Sampling', linewidth = 0.5, edgecolor='black', ax=ax)
stations_used[0:15].plot(x="Lon", y="Lat", marker= '^', kind="scatter", s = 200, color = 'yellow', linewidth = 0.5, edgecolor='black', label = 'USGS/EPA Gauges', ax=ax)

# add grid
ax.set_xticks([-112, -111])
ax.grid(visible=True, alpha=0.5)
ax.get_legend().remove()
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
plt.show()

################ Figure 4: All data points vs GMWL ################

fig, ax = plt.subplots(figsize=(12,12))
x = np.linspace(-25,0,100)
y = 8*x + 10
plt.plot(x, y, color = 'black', label='GMWL')
bl.loc[(bl['Subgroup'] == "Out")].plot(x="d18O", y="dD", kind='scatter', color = '#3D7BBA', marker='P', s=150, linewidth = 0.5, edgecolor='black', ax=ax, label = 'Out')
bl.loc[(bl['Subgroup'] == "Downstream")].plot(x="d18O", y="dD", kind='scatter', color = '#F8C537', marker='v', s=100, linewidth = 0.5, edgecolor='black', ax=ax, label = 'Downstream')
bl.loc[(bl['Subgroup'] == "Upstream")].plot(x="d18O", y="dD", kind='scatter', color = '#FF7F50', marker='^', s=100, linewidth = 0.5, edgecolor='black', ax=ax, label = 'Upstream')
bl.loc[(bl['Subgroup'] == "Around")].plot(x="d18O", y="dD", kind='scatter', color = '#00C176', marker='D', s=70, linewidth = 0.5, edgecolor='black', ax=ax, label = 'Data around lake')
ax.set_xlabel('δ$^1$$^8$O (‰)')
ax.set_ylabel('δ$^2$H (‰)')
plt.show()

################ Figure 5: Plot changing x and humidities against GMWL ################

## Inputs:
    
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

# Input X values derived from d18O and dD mass balance calculations 1

x_O = 0.495087888113281 #d18O
x_D = 0.368869387932284 #dD

## Generate array of X and humidity values

x_arr = np.arange(0,1.1,0.1)
humidity = [0.0,0.2,0.4,0.6,0.76,0.95]

## Plot

fig, ax = plt.subplots(figsize=(9,10))

x = np.linspace(-20,5,100)
y = 8*x + 10 # GMWL
ax.plot(x, y, color = 'black')

### Extract current lake composition and theoretical maximum of lake enrichment:
### Function mass_balance_ssx output: [dX_S (lake isotope), limit (theoretical maximum enrichment)]
    
LS_Ox = lbf.mass_balance_ssx(hum, ep_k_O, ep_O, alpha_O, atm_O, influx_O, x_O) # should be similar to lake_O
LS_Dx = lbf.mass_balance_ssx(hum, ep_k_D, ep_D, alpha_D, atm_D, influx_D, x_D) # should be similar to lakd_D

for humx in humidity:
    ep_k2_O = lbf.kinetic_en_d18O(humx)
    ep_k2_D = lbf.kinetic_en_dD(humx)
    
    dX_LS_O = []
    dX_LS_D = []
    for xx in x_arr:
        LS_O, l = lbf.mass_balance_ssx(humx, ep_k2_O, ep_O, alpha_O, atm_O, influx_O, xx)
        dX_LS_O.append(LS_O)
        LS_D, l = lbf.mass_balance_ssx(humx, ep_k2_D, ep_D, alpha_D, atm_D, influx_D, xx)
        dX_LS_D.append(LS_D)
    ax.plot(dX_LS_O,dX_LS_D, color='grey')
#ax.scatter(-9.147563681404433, -83.53181856351087, marker="D", color = 'black', s=150, edgecolor='black', zorder=10, label = "Back-calculated lake composition")
ax.scatter(LS_Ox[0],  LS_Dx[0], marker="D", color = 'red', s=150, edgecolor='black', zorder=10, label = "Current lake isotopic composition")
ax.scatter(LS_Ox[1],  LS_Dx[1], marker="P", color = 'black', s=200, edgecolor='black', zorder=10, label = "Theoretical maximum enrichment")

bl.loc[(bl['Subgroup'] == "Around")].plot(x="d18O", y="dD", kind='scatter', edgecolor='black', color = 'white', label = 'Data around lake', ax=ax)
plt.legend(bbox_to_anchor=(1, 1.0), fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.xlabel("δ$^1$$^8$O (‰)", fontsize=15)
plt.ylabel("δ$^2$H (‰)", fontsize=15)
plt.show()

################ Figure 7: Mean isotopic composition of different components ################

fig, ax = plt.subplots(figsize=(16,16))
x = np.linspace(-25,0,100)
y = 8*x + 10
plt.plot(x, y, color = 'black', label='GMWL', zorder=0)
bl.loc[(bl['Subgroup'] == "Around")].plot(x="d18O", y="dD", kind='scatter', color = 'white', marker='o', s=70, linewidth = 0.5, edgecolor='black', alpha=0.5, ax=ax, label = 'Data around lake')
ax.scatter(-28.22, -179.88, color='black', marker = '*', linewidth = 1, edgecolor='black', s=700, label='Evaporate (δE)')
ax.scatter(-22.08, -163.81, color='white', marker = '*', linewidth = 1, edgecolor='black', s=700, label='Atmosphere (δA)')
ax.scatter(-14.6, -105.7, color='white', marker = '^', linewidth = 1, edgecolor='black', s=700, label='Volume-weighted precipitation (δP)')
ax.scatter(-11.7, -84.02, color='gray', marker = '^', linewidth = 1, edgecolor='black', s=700, label='Evaporation-flux weighted precipitation (δP)')
ax.scatter(-16.22, -122.15, color='black', marker = '^', linewidth = 1, edgecolor='black', s=700, label='Total inflow (δI)')
ax.scatter(-8.76, -86.44, color='black', marker = 'D', linewidth = 1, edgecolor='black', s=400, label='Steady-state lake (δS)')
ax.plot([-28.22, -8.76, -16.22], [-179.88, -86.44, -122.15], ls = '--', zorder=0)

plt.legend(labelspacing = 1)
ax.set_xlabel('δ$^1$$^8$O (‰)')
ax.set_ylabel('δ$^2$H (‰)')
plt.show()
