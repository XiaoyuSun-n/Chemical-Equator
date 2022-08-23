# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:59:21 2021

@author: 56907
"""

import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.util import add_cyclic_point
from palettable.colorbrewer.diverging import RdBu_11_r
import scipy.interpolate
import datetime
import heapq

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def days(str1,str2):
#To calculate the days between to DATE (in the formate %Y-%m-%d)
#date1 needs to be after date2
    date1=datetime.datetime.strptime(str1[0:10],'%Y-%m-%d')
    date2=datetime.datetime.strptime(str2[0:10],'%Y-%m-%d')
    num=(date1-date2).days
    return num

###
#The grid (latitude and longitude) in a fine resolution 0.5 x 0.625 degrees of the simulation cases
#
latitute=np.load('D:\\Geos_output\\py\\data\\nest\\latitute_nest.npy')
longitute=np.load('D:\\Geos_output\\py\\data\\nest\\longitute_nest.npy')

latitute_c=np.load('F:\\ITCZ\\nh\\latitude.npy')
longitute_c=np.load('F:\\ITCZ\\nh\\longitude.npy')

#
#The simulation results from Geos-chem. Can be downloaded from the link. Details please find in README 
#Ptracer_mn[day,level,latitude,longitude]
#Note: The calculation needs the simulation with a fine resolution, 0.5 x 0.625 degrees.
#      So data from '201912.nc4' is in a fine resolution.
#
AA      = xr.open_dataset("F:\\passivetracer\\w_h\\sh\\201912.nc4")
Ptracer_mn=AA['SpeciesConc_PassiveTracer'].data

###
#Trend component from the 'intercept_finding.py'. It needed to be calculated and saved after running the script.
#trend[level,day] 
#Note: the trend is five year decompose results, so it should be treated properly by the function 'days'. 
#
trend=np.load('F:\\ITCZ\\vertical\\sh\\mid_vl.npy')

#
#'2019-12-01' is the month which the Chemical Equator (CE) is calculated. 
#This infomation should be changed in each month of interest.
# 
#e.g. '2019-12-01' means that it calculated the CE in December 2019,
#and this should be consistent with the data loaded before ('201912.nc4')
#

trend_days=days('2019-12-01','2015-01-01')

ce=np.zeros((Ptracer_mn.shape[0],trend.shape[0],Ptracer_mn.shape[3]))
for k in range(Ptracer_mn.shape[0]):
    for i in range(trend.shape[0]):
        for j in range(Ptracer_mn.shape[3]):
            ce_idx=find_nearest(Ptracer_mn[k,i,:,j],trend[i,k+trend_days])
            ce[k,i,j]=latitute[ce_idx]
        

###
diff=np.empty(ce.shape,dtype=float)

for k in range (ce.shape[0]):
    for i in range (ce.shape[1]):
        for j in range (ce.shape[2]-1):
            diff[k,i,j]=ce[k,i,j+1]-ce[k,i,j]


for i in range(diff.shape[0]):
    for j in range(diff.shape[1]):
        dif_std=np.std(diff[i,j,:])
        dif_ave=np.mean(diff[i,j,:])
        dif_out=np.where(((3*dif_std+dif_ave)<diff[i,j,:]) | ((-3*dif_std+dif_ave)>diff[i,j,:]) )[0]
        dif_in=np.arange(0,ce.shape[2],1)
        dif_in[dif_out]=-9999
        print('i=',i)
        for k in range(dif_out.shape[0]):
            print('j=',k)
            delta=abs(dif_out[k]-dif_in)
            x_f_min=np.where(delta==heapq.nsmallest(2,delta)[0])[0][0]
            print(x_f_min)
            ce[i,j,dif_out[k]]=ce[i,j,x_f_min]

###            
#a preliminary plot for visializing the horizontal result of CE in a level of interest       
plt.figure()
projection = ccrs.PlateCarree()        
fig, ax = plt.subplots(figsize=(16, 9), subplot_kw=dict(projection=projection))
#Ptracer_mn averaged in month, the second axis is the level which can be changed of interest
cs = ax.contourf(longitute[:-1], latitute[:-1], np.mean(Ptracer_mn[:,30,:,:]*1e3,axis=0), cmap=RdBu_11_r.mpl_colormap,levels=15)
pitcz = plt.plot(longitute[:-1],np.mean(ce[0:30,0,:],axis=0),color='k',linewidth=2)
ax.coastlines(linewidth=1.5) 
ax.set_xticks(np.arange(-180, 181, 60), crs=projection)
ax.set_yticks(np.arange(-30, 31, 30), crs=projection)
# set ticklabels
lon_formatter = LongitudeFormatter(number_format='.0f',degree_symbol='',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.0f',degree_symbol='')
cb = fig.colorbar(cs, orientation='horizontal',shrink=0.73, pad=0.05)   
fig.set_size_inches(16,9)

#a preliminary plot for visializing the vertical result of CE in a longitude range of interest  
plt.figure()
projection = ccrs.PlateCarree()
#Ptracer_mn averaged in longitude range, the second axis is the level which can be changed of interest        
cs = plt.contourf(latitute[:-1], AA.lev[0:35],np.mean(Ptracer_mn[20,0:35,:,480:575]*1e3,axis=2), cmap=RdBu_11_r.mpl_colormap,levels=30)
pitcz = plt.plot(np.mean(ce[20,:,480:575],axis=1),AA.lev[0:35],color='k',linewidth=2)
plt.gca().invert_yaxis()
cb = fig.colorbar(cs, orientation='horizontal',shrink=0.73, pad=0.05)   



















