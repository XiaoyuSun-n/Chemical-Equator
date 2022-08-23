# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 13:01:04 2021

@author: 56907
"""


import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os,time
from datetime import datetime
from matplotlib.dates import date2num
import tables
from scipy import interpolate
import seaborn as sns
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
from palettable.colorbrewer.diverging import RdBu_11_r 
import netCDF4 as nc
from numpy import polyfit, poly1d
from statsmodels.tsa.seasonal import seasonal_decompose
import time

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#
#To load the simulation results.
#'daily_20xx.nc4' contains the simulation results of the passive tracer in the year 20xx.
#This is an example of the northern hemisphere cases
#
    
nh2015       = xr.open_dataset("F:\\passivetracer\\sh\\daily_2015.nc4")
nh2016       = xr.open_dataset("F:\\passivetracer\\sh\\daily_2016.nc4")
nh2017       = xr.open_dataset("F:\\passivetracer\\sh\\daily_2017.nc4")
nh2018       = xr.open_dataset("F:\\passivetracer\\sh\\daily_2018.nc4")
nh2019       = xr.open_dataset("F:\\passivetracer\\sh\\daily_2019.nc4")

##

start=time.time()
lat_range = np.where( (nh2015.lat <=30) & (nh2015.lat>=-30))[0]
lev_range = np.where( (nh2015.lev <=1) & (nh2015.lev>=0.1))[0]
nhyear_tropical = xr.combine_nested([nh2015['SpeciesConc_PassiveTracer'].isel(lat=lat_range,lev=lev_range),
                            nh2016['SpeciesConc_PassiveTracer'].isel(lat=lat_range,lev=lev_range),
                            nh2017['SpeciesConc_PassiveTracer'].isel(lat=lat_range,lev=lev_range),
                            nh2018['SpeciesConc_PassiveTracer'].isel(lat=lat_range,lev=lev_range),
                            nh2019['SpeciesConc_PassiveTracer'].isel(lat=lat_range,lev=lev_range)],
                            concat_dim=['time'])
end=time.time()
print (str(end-start))

##

season_tracer = np.zeros((len(nhyear_tropical.lev),len(nhyear_tropical.lat),len(nhyear_tropical.lon),len(nhyear_tropical.time)))
residual_tracer = np.zeros((len(nhyear_tropical.lev),len(nhyear_tropical.lat),len(nhyear_tropical.lon),len(nhyear_tropical.time)))
trend_tracer = np.zeros((len(nhyear_tropical.lev),len(nhyear_tropical.lat),len(nhyear_tropical.lon),len(nhyear_tropical.time)))

start=time.time()
for i in range (len(nhyear_tropical.lat)):
    for j in range (len(nhyear_tropical.lon)):
        for k in range (len(nhyear_tropical.lev)):          
            analysis=nhyear_tropical[:,k,i,j]
            decompose_result_mult = seasonal_decompose(analysis,period=365, model="additive")
            trend = decompose_result_mult.trend
            seasonal = decompose_result_mult.seasonal
            residual = decompose_result_mult.resid
            season_tracer[k,i,j,:]   = seasonal
            trend_tracer[k,i,j,:] = trend
            residual_tracer[k,i,j,:]=residual
    print(i)
end=time.time()
print (str(end-start))

start=time.time()
x1=np.arange(0,trend_tracer[0,0,0,~np.isnan(trend_tracer[0,0,0,:])].shape[0],1)+182
x2=np.arange(0,trend_tracer.shape[3],1)
trend2 = np.zeros((trend_tracer.shape[0],trend_tracer.shape[1],trend_tracer.shape[2],trend_tracer.shape[3]))
for i in range(trend_tracer.shape[0]):
    for j in range (trend_tracer.shape[1]):
        for k in range (trend_tracer.shape[2]):
            z1 = np.polyfit(x1,trend_tracer[i,j,k,~np.isnan(trend_tracer[i,j,k,:])], 1)
            trend2[i,j,k,:]=z1[0]*x2+z1[1]
    print(i)
end=time.time()
print (str(end-start))

start=time.time()
mid_vl=np.zeros((trend2.shape[0],trend2.shape[3]))
for i in range(mid_vl.shape[0]):
    for j in range(mid_vl.shape[1]):
        mid_vl[i,j]=np.mean(trend2[i,:,:,j])
    print(i)
end=time.time()
print (str(end-start))   

## This variable 'mid_vl' needs to be saved, which is the trend component used in 'chemical_equator_multi_level.py'





    
    
    
    
    
    



