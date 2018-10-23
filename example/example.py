# -*- coding: utf-8 -*-

import numpy as np
from epistemic_uncertainties import uncertainties as uncer
from epistemic_uncertainties import changeFault as chf

faultdir = './'
faultfile = 'fault_patches_14_11.rectangles'
nstrike = 14
utmzone = 32
dip_range = [0,1,2]

# Initialize the class changeFault
faultvar = chf(faultdir,faultfile,nstrike,utmzone,inputCoordinates='lonlatz', nbrcorners=4)

# Write new fault geometries with a dip deeper of 1 and 2 degrees
faultvar.varyDip(dip_range)
faultvar.varyStrike([1,2],'center')
 
Np = 308  # number of parameters
GFs = ['G_dip_'+str(d)+'.txt' for d in dip_range]  # list of the Green's functions
sigma_dip = 5  # standard deviation on the dip parameter
mprior = np.ones((Np))

uncert = uncer(faultdir, Np, multifault=False, export='./')
uncert.calcCpDip(GFs,dip_range,sigma_dip,mprior)
