#!/usr/bin/env python-atmos
# conda execute
# env:
#  - python >=3
#  - numpy

import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
import os.path
import sys


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpl_toolkits
from  mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA



''' Constants '''
a   = 6371e3      # m
g   = 9.81      # m/s
rho = 1.225     # kg m^-3
Rd  = 287.05    # J Kg^-1 K^-1
Rv  = 461.51    # J Kg^-1 K^-1
Cpd = 1005      # J Kg^-1 K^-1
Cvd = 718.0
Cpv = 1880      # J Kg^-1 K^-1
Lv = 2510.4e3   # J Kg^-1
Tref = 273.15
R_rat = Rv / Rd - 1.0

CO2 = '300'
S0  = '1151'
var = ['temp', 'ps','sphum','ucomp','vcomp']
mtime = [ 248, 224, 248, 240, 248, 240, 248, 248, 240, 248, 240, 248 ]
print(mtime)
print(mtime[2])

nmonth = 12
dirname = "/scratch/hkb2/HIRAM_AQUA/aqua/aqua_som20_qflux0_"+CO2+"co2"\
          "_delsol1.0_cosz2.0_S"+S0+"_alb8_zspecdaily_mean_heat8.e07/output/" 


''' dp '''
a32 = np.array([100.00000,     400.00000,     818.60211,\
	1378.88653,    2091.79519,    2983.64084,\
	4121.78960,    5579.22148,    7419.79300, \
	9704.82578,   12496.33710,   15855.26306, \
	19839.62499,   24502.73262,   28177.10152, \
	29525.28447,   29016.34358,   27131.32792, \
	24406.11225,   21326.04907,   18221.18357, \
	15275.14642,   12581.67796,   10181.42843, \
	8081.89816,    6270.86956,    4725.35001, \
	3417.39199,    2317.75459,    1398.09473, \
	632.49506,       0.00000,       0.00000] )
b32 = np.array([0.00000,       0.00000,       0.00000, \
	0.00000,       0.00000,       0.00000, \
	0.00000,       0.00000,       0.00000, \
	0.00000,       0.00000,       0.00000, \
	0.00000,       0.00000,       0.01711, \
	0.06479,       0.13730,       0.22693, \
	0.32416,       0.42058,       0.51105, \
	0.59325,       0.66628,       0.73011, \
	0.78516,       0.83217,       0.87197, \
	0.90546,       0.93349,       0.95685, \
	0.97624,       0.99223,       1.00000] )
phalf = np.zeros([33,360,720])
dp   = np.zeros([32,360,720])
phi  = np.zeros([32,360,720]) 
DSE_trans_tot = np.zeros([nmonth,360])
DSE_trans_eul = np.zeros([nmonth,360])
SH_trans_tot = np.zeros([nmonth,360])
SH_trans_eul = np.zeros([nmonth,360])
LH_trans_tot = np.zeros([nmonth,360])
LH_trans_eul = np.zeros([nmonth,360])
time = -1
for month in range(nmonth):
    print('month =', month)
    time = time+1
    for time in range(time,time+mtime[month]): 
        print('time =', time)
        for j in var:
            filename ="postproc1995.atmos_8xdaily_aqua_som20_qflux0_"+CO2+"co2_delsol1.0_"\
                      "cosz2.0_S"+S0+"_alb8_zspecdaily_mean_heat8.e07_"+j+"1995.nc"
            f = Dataset(os.path.join(dirname+filename))
            if j is 'ps':
                vars()[j] = np.array(f.variables[j][time,:,:])
                for k in range(33):
                    phalf[k,:,:] = a32[k] + b32[k] * ps
                    dp = phalf[1:33,:,:] - phalf[0:32,:,:]
            else:
                vars()[j] = np.array(f.variables[j][time,:,:,:])
        lat   =  np.array(f.variables['grid_yt'])
        vars()["vcorr"] = np.sum((vcomp * dp), axis=(0,2)) / np.sum( dp, axis=(0,2) )
        vars()["v"] = vcomp + vcorr[np.newaxis,:,np.newaxis]
        vars()["tv"]  = temp * (1+ R_rat * sphum)
        phi[31,:,:] = 0.5 * dp[31,:,:] * Rd * tv[31,:,:] /phalf[32,:,:] 
        for k in range(30,0,-1):
              phi[k,:,:] = phi[k+1,:,:] + 0.5 * dp[k,:,:] * Rd * tv[k,:,:] /phalf[k,:,:]\
              + 0.5 * dp[k+1,:,:] * Rd * tv[k+1,:,:] /phalf[k+1,:,:]
        
        DSE_trans_tot[month,:] = DSE_trans_tot[month,:] + 2 * np.pi * a * np.cos(np.deg2rad(lat)) \
                                 * np.sum(np.mean(v*( Cpd * (temp-300) + phi) * dp,2),0) / g
        DSE_trans_eul[month,:] = DSE_trans_eul[month,:] + 2 * np.pi * a * np.cos(np.deg2rad(lat)) \
                                 * np.sum(np.mean(v*dp,2) * Cpd * np.mean((temp-300),2) + np.mean(phi,2), 0) / g
        SH_trans_tot[month,:]  = SH_trans_tot[month,:] + 2 * np.pi * a * np.cos(np.deg2rad(lat))\
                                 * np.sum(np.mean(v* Cpd * (temp-300) * dp,2),0) /g
        SH_trans_tot[month,:]  = SH_trans_tot[month,:] + 2 * np.pi * a * np.cos(np.deg2rad(lat))\
                                 * np.sum(np.mean(v*dp,2) * Cpd * np.mean((temp-300), 2), 0) /g

        LH_trans_tot[month,:]  = LH_trans_tot[month,:] + 2 * np.pi * a * np.cos(np.deg2rad(lat))\
                                 * np.sum(np.mean(v* Lv * sphum * dp,2),0) /g
        LH_trans_eul[month,:]  = LH_trans_eul[month,:] + 2 * np.pi * a * np.cos(np.deg2rad(lat))\
                                 * np.sum(np.mean(v*dp,2) * Lv * np.mean(sphum, 2), 0) /g
        


data_out= np.vstack((DSE_trans_tot, DSE_trans_eul, SH_trans_tot, SH_trans_eul, LH_trans_tot, LH_trans_eul )).T
filename_out = "Transport_test.dat"
np.savetxt(filename_out,data_out, delimiter ="\t")


