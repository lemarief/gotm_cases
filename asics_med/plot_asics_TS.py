import xarray as xr
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import pandas as pd

class asics_case:
    #============================================================================
    #==
    #============================================================================    
    def __init__(self,data_params):
        default_params = { "filename": "dataset.nc", "title": "toto"}
        params = default_params.copy()
        params.update(data_params)
        self.config = params
        self.title  = self.config["title"]
    
    def get_data_gotm(self):
       ds = xr.open_dataset(self.config["filename"])
       self.time   = xr.decode_cf(ds).time
       self.temp   = ds["temp_p"].squeeze()   # (time, z)
       self.salt   = ds["salt_p"].squeeze()
       self.v      = ds["v"].squeeze()
       self.tke    = ds["tke"].squeeze()
       depth       = ds["z"].squeeze()
       depthi      = ds["zi"].squeeze()
       self.z      = depth.mean(dim="time")
       self.zi     = depthi.mean(dim="time")

    def get_data_obs(self):
       ds = xr.open_dataset(self.config['filename']) 
       ds = ds.assign_coords(time_counter = pd.to_datetime("2013-01-01") + pd.to_timedelta(ds.time_counter, unit="D"))
       ds_sel = ds.sel(time_counter=slice("2013-01-15", "2013-02-10")) 
       self.temp = ds_sel["votemper"]
       self.salt = ds_sel["vosaline"]
       self.time = ds_sel["time_counter"] 
       self.z    = -1. * ds["deptht"].squeeze()

    def plot_temp(self,ax,tmin,tmax,colorMap,Tlevels):
       pcm  = ax.pcolormesh(self.time, self.z, self.temp.T,shading="auto",cmap=colorMap,vmin=tmin,vmax=tmax)  
       cs   = ax.contour(self.time, self.z, self.temp.T, colors='k', levels=Tlevels, linewidths=1.5)
       ax.clabel(cs, inline=True, fontsize=10, fmt = '%1.2f')                        
       ax.set_ylabel("Depth (m)")
       ax.set_title(self.title)
       plt.colorbar(pcm, ax=ax, label="°C")

    def plot_salt(self,ax,smin,smax,colorMap,Slevels):
       pcm  = ax.pcolormesh(self.time, self.z, self.salt.T,shading="auto",cmap=colorMap,vmin=smin,vmax=smax)  
       cs   = ax.contour(self.time, self.z, self.salt.T, colors='k', levels=Slevels, linewidths=1.5)
       ax.clabel(cs, inline=True, fontsize=10, fmt = '%1.2f')                        
       ax.set_ylabel("Depth (m)")
       ax.set_title(self.title)
       plt.colorbar(pcm, ax=ax, label="psu")

    def plot_yvel(self,ax,smin,smax,colorMap,Slevels):
       pcm  = ax.pcolormesh(self.time, self.z, self.v.T,shading="auto",cmap=colorMap,vmin=smin,vmax=smax)  
       cs   = ax.contour(self.time, self.z, self.v.T, colors='k', levels=Slevels, linewidths=1.5)
       ax.clabel(cs, inline=True, fontsize=10, fmt = '%1.2f')                        
       ax.set_ylabel("Depth (m)")
       ax.set_title(self.title)
       plt.colorbar(pcm, ax=ax, label="m/s")

params = {'filename': 'obs_data/AsicsMed_data.nc','title': 'Observational data'}
scm_obs = asics_case(params)
scm_obs.get_data_obs()


runs = [ {'filename': 'Asics_Med_keps.nc','title': r'$k-\epsilon$'},
         {'filename': 'Asics_Med_edmf.nc','title': 'EDMF'},
         {'filename': 'Asics_Med_edmf_energy.nc','title': 'EDMF + MF on dynamics + Energy'}
]

nrun = len(runs)
scm  = [0]*len(runs)
for i, run_params in enumerate(runs):
   scm[i] = asics_case(run_params)
   scm[i].get_data_gotm() 

#==
tmin=12.875; tmax=13.205; colorMapT='twilight'; titleT='Temperature [Celsius]'
smin=38.45 ; smax=38.55 ; colorMapS='viridis' ; titleS='Salinity [Celsius]'
hmax=2300. ; Tlevels=[12.90, 13, 13.10]; Slevels=[38.49, 38.50, 38.52]
#== 

# --------------------------------------------------
# Création des figures
# --------------------------------------------------
fig, ax = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
#
scm_obs.plot_temp(ax[0],tmin,tmax,colorMapT,Tlevels)
scm[0].plot_temp(ax[1],tmin,tmax,colorMapT,Tlevels)
scm[1].plot_temp(ax[2],tmin,tmax,colorMapT,Tlevels)
scm[2].plot_temp(ax[3],tmin,tmax,colorMapT,Tlevels)

fig, ax = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
#
scm_obs.plot_salt(ax[0],smin,smax,colorMapS,Slevels)
scm[0].plot_salt(ax[1],smin,smax,colorMapS,Slevels)
scm[1].plot_salt(ax[2],smin,smax,colorMapS,Slevels)
scm[2].plot_salt(ax[3],smin,smax,colorMapS,Slevels)

fig, ax = plt.subplots(3, 1, figsize=(12, 7.5), sharex=True)
#
vmin = -0.025
vmax =  0.025
Vlevels = [-0.2, 0., 0.02]
colorMapV='RdBu'
scm[0].plot_yvel(ax[0],vmin,vmax,colorMapV,Vlevels)
scm[1].plot_yvel(ax[1],vmin,vmax,colorMapV,Vlevels)
scm[2].plot_yvel(ax[2],vmin,vmax,colorMapV,Vlevels)

plt.tight_layout()
plt.show()