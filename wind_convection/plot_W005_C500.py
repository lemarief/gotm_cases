#!/usr/bin/env python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
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
        nc = Dataset(self.config["filename"] , mode='r')
        self.temp = nc.variables["temp"][:]
        self.z    = nc.variables["z"][:]
        self.tke  = nc.variables["tke_res"][:]+nc.variables["tke_sbg"][:]
        self.wt   = nc.variables["wt_res"][:]+nc.variables["wt_sbg"][:]
        self.wtke = nc.variables["wtke_res"][:]+nc.variables["wtke_sbg"][:]
        self.zinv = - self.z[np.argmin(self.wt)]
        if self.config["momentum_diag"]:
            self.u  = nc.variables["u"][:]
            self.wu = nc.variables["wu_res"][:]+nc.variables["wu_sbg"][:]
        nc.close()
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
                + 0.5*nc.variables["w_p"][-1, 1:, 0, 0]*nc.variables["w_p"][-1, 1:, 0, 0] ) 

        if self.config["momentum_diag"]:
            self.u  = nc.variables["u"][-1, :, 0, 0]
            akm     = nc.variables["num"][-1, :, 0, 0]
            self.wu = - akm[1:-1]*( self.u[1:]-self.u[:-1] ) / ( self.z[1:]-self.z[:-1] )
            if self.config["mass_flux"]:
               self.wu = self.wu - nc.variables["fmass"][-1, 1:-1, 0, 0]*( nc.variables["u_p"][-1, 1:-1, 0, 0] - self.u[1:] )  
            if self.config["mass_flux_energy"]:
               self.wtke = self.wtke - nc.variables["fmass"][-1, 1:, 0, 0]*(
                            0.5*(nc.variables["u_p"][-1, 1:, 0, 0]-self.u[:]   )**2 ) 
        nc.close()
    #============================================================================
    #==
    #============================================================================
runs = [ {'filename': 'W005_C500_keps.nc', 'momentum_diag': True},
         {'filename': 'W005_C500_tke.nc', 'momentum_diag': True},
         {'filename': 'W005_C500_edmf.nc', 'momentum_diag': True, 'mass_flux': True},
         {'filename': 'W005_C500_edmf_energy.nc', 'momentum_diag': True, 'mass_flux': True, "mass_flux_energy": True},
]

nrun = len(runs)
scm  = [0]*len(runs)
for i, run_params in enumerate(runs):
   scm[i] = convection(run_params)
   scm[i].get_data_gotm() 

params = {'filename': 'les_data/W005_C500_LES_t72h.nc', 'momentum_diag': True, 'mass_flux': False}
scm_les = convection(params)
scm_les.get_data_les()

vscale = 1./scm_les.zinv
zmax   = -1.2

colors = ['tab:orange','tab:blue','tab:green','tab:green']
styles = ['solid', 'solid', 'dashed','solid']
linew  = [2,2,3,3]
label  = [r"$k-\epsilon$",r"TKE-1eq",r"EDMF",r"EDMF-Energy"]

fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False,sharey=True, figsize=(12, 8) )

# ROW1
ax = axes.flat[0]
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

ax = axes.flat[1]
ax.set_title(r"$\overline{u}$",fontsize=18)
ax.plot(scm_les.u, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].u, vscale*scm[i].z, linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
xmin = 0.; xmax = 0.075
ax.set_xlim(xmin,xmax)
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"${\rm m}\;{\rm s}^{-1}$",fontsize=18)
ax.set_ylabel(r"$z/h$",fontsize=14)

ax = axes.flat[2]
ax.set_title(r"$k$",fontsize=18)
ax.plot(scm_les.tke, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].tke[:-1], vscale*scm[i].zi[:-1], linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$m^{2}\;s^{-2}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# ROW2
ax = axes.flat[3]
ax.set_title(r"$\overline{w^\prime \theta^\prime}$",fontsize=18)
ax.plot(scm_les.wt, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].wt, vscale*scm[i].zi[1:-1], linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$^{\circ}{\rm C}\;m\;s^{-1}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


ax = axes.flat[4]
ax.set_title(r"$\overline{w^\prime u^\prime}$",fontsize=18)
ax.plot(scm_les.wu, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].wu, vscale*scm[i].zi[1:-1], linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$m^{2}\;s^{-2}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))


ax = axes.flat[5]
ax.set_title(r"$w^\prime k$",fontsize=18)
ax.plot(scm_les.wtke, vscale*scm_les.z, 'ko', linewidth=2.0, label='LES')
for i in range(nrun):
    ax.plot(scm[i].wtke[:], vscale*scm[i].z[:], linewidth=linew[i], color = colors[i], linestyle=styles[i], label = label[i])
#==
ax.set_ylim(zmax, 0)
ax.set_xlabel(r"$m^{3}\;s^{-3}$",fontsize=18)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

fig.tight_layout()
plt.show()
