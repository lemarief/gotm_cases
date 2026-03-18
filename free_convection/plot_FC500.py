#!/usr/bin/env python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#
plt.rcParams['text.usetex'] = True

class convection:
    #============================================================================
    #==
    #============================================================================    
    def __init__(self,data_params):
        default_params = { "filename": "dataset.nc",  "momentum_diag": False, "mass_flux": False, "mass_flux_energy": False }
        params = default_params.copy()
        params.update(data_params)
        self.config = params
    #============================================================================
    #==
    #============================================================================
    def get_data_les(self):
        df = pd.read_csv(self.config["filename"], sep="\t") 
        self.z    = df["z"].values
        self.tke  = df["tke_res"].values + df["tke_sbg"].values
        self.temp = df["temp"].values
        self.wt   = df["wt_res"].values+df["wt_sbg"].values
        self.wtke = df["wtke_res"].values+df["wtke_sbg"].values
        self.zinv = - self.z[np.argmin(self.wt)]
        if self.config["momentum_diag"]:
            self.u  = df["u"].values
            self.wu = df["wu_res"].values+df["wu_sbg"].values
    #============================================================================
    #==
    #============================================================================
    def get_data_gotm(self):
        nc = Dataset(self.config["filename"] , mode='r')
        if "temp_p" in nc.variables:
           self.temp = nc.variables["temp_p"][-1, :, 0, 0]
        else:
           self.temp = nc.variables["temp"][-1, :, 0, 0]
        self.z    = nc.variables["z"][-1, :, 0, 0]
        self.zi   = nc.variables["zi"][-1, :, 0, 0]
        self.tke  = nc.variables["tke"][-1, :, 0, 0] 
        akh       = nc.variables["nuh"][-1, :, 0, 0]
        self.wt   = - akh[1:-1]*( self.temp[1:]-self.temp[:-1] ) / ( self.z[1:]-self.z[:-1] )   

        if self.config["mass_flux"]:
            self.wt = self.wt - nc.variables["fmass"][-1, 1:-1, 0, 0]*( nc.variables["t_p"][-1, 1:-1, 0, 0] - self.temp[1:] )  
        ake       = 0.5*(nc.variables["nuh"][-1, :-1, 0, 0]+nc.variables["nuh"][-1, 1:, 0, 0])
        self.wtke = - ake *( self.tke[1:]-self.tke[:-1] ) / ( self.zi[1:]-self.zi[:-1] ) 

        if self.config["mass_flux"] and self.config["mass_flux_energy"]:
            self.wtke = self.wtke - nc.variables["fmass"][-1, 1:, 0, 0]*( 
                ( nc.variables["tke_p"][-1, 1:, 0, 0]-nc.variables["tke"][-1, :-1, 0, 0] ) 
                + 0.5*nc.variables["w_p"][-1, 1:, 0, 0]*nc.variables["w_p"][-1, 1:, 0, 0]  ) 

        if self.config["momentum_diag"]:
            self.u  = nc.variables["u"][-1, :, 0, 0]
            akm     = nc.variables["num"][-1, :, 0, 0]
            self.wu = - akm[1:-1]*( self.u[1:]-self.u[:-1] ) / ( self.z[1:]-self.z[:-1] )

        nc.close()
    #============================================================================
    #==
    #============================================================================
runs = [ {'filename': 'FC500_keps.nc', 'momentum_diag': False, 'mass_flux': False, "mass_flux_energy": False},
         {'filename': 'FC500_tke.nc', 'momentum_diag': False, 'mass_flux': False, "mass_flux_energy": False},
         {'filename': 'FC500_edmf.nc', 'momentum_diag': False, 'mass_flux': True, "mass_flux_energy": False},
         {'filename': 'FC500_edmf_energy.nc', 'momentum_diag': False, 'mass_flux': True, "mass_flux_energy": True},
]

nrun = len(runs)
scm  = [0]*len(runs)
for i, run_params in enumerate(runs):
   scm[i] = convection(run_params)
   scm[i].get_data_gotm() 

params = {'filename': 'les_data/FC500_LES_t72h.dat', 'momentum_diag': False, 'mass_flux': False}
scm_les = convection(params)
scm_les.get_data_les()

vscale = 1./scm_les.zinv
zmax   = -1.2

colors = ['tab:orange','tab:blue','tab:green','tab:green']
styles = ['solid', 'solid', 'dashed','solid']
linew  = [2,2,3,3]
label  = [r"$k-\epsilon$",r"TKE-1eq",r"EDMF",r"EDMF-Energy"]

fig, axes = plt.subplots(nrows=1, ncols=4, sharex=False,sharey=True, figsize=(12, 5) )


ax = axes[0]
ax.set_title(r"$\overline{\theta}$",fontsize=18)
ax.plot(scm_les.temp, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].temp, vscale*scm[i].z, linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
xmin = 1.65; xmax = 1.80
ax.set_xlim(xmin,xmax)
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$^{\circ}{\rm C}$",fontsize=18)
ax.set_ylabel(r"$z/h$",fontsize=14)
ax.legend(loc=6,fontsize=10)


ax = axes[1]
ax.set_title(r"$\overline{w^\prime \theta^\prime}$",fontsize=18)
ax.plot(scm_les.wt, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].wt, vscale*scm[i].zi[1:-1], linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$^{\circ}{\rm C}\;m\;s^{-1}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


ax = axes[2]
ax.set_title(r"$k$",fontsize=18)
ax.plot(scm_les.tke, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].tke[:-1], vscale*scm[i].zi[:-1], linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$m^{2}\;s^{-2}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

ax = axes[3]
ax.set_title(r"$w^\prime k$",fontsize=18)
ax.plot(scm_les.wtke, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].wtke[:], vscale*scm[i].z[:], linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$m^{2}\;s^{-2}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

fig.tight_layout()
plt.show()
