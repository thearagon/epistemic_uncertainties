## epistemic_uncertainties
Account for uncertainties in the assumed fault geometry. The purpose of this module is to calculate covariance matrices of the predictions derived from the uncertainties in the assumed fault geometry.

This methodology is used in [Ragon et al. (2018)]

  [Ragon et al. (2018)]: https://academic.oup.com/gji/advance-article/doi/10.1093/gji/ggy187/4996353?guestAccessKey=86188bcd-f078-4700-9762-27d888995596
 
To account for uncertainties in the Earth elastic properties, please refer to [Duputel et al. (2014)]

  [Duputel et al. (2014)]: https://academic.oup.com/gji/article/197/1/464/686313

This module consists in two classes:
- changeFault: This class helps modifying the assumed fault geometry to calculate various Green's functions (the tool to calculate Green's functions is not furnished here).
- uncertainties: This class allows us to calculate the uncertainties in the predictions.

In the following, you will find:
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

## Installation
No installation needed. Just download the file epistemic_uncertainties.py 

## Usage

To compute various fault geometries to calculate Green's functions:

1) Write your fault geometry in a file with a GMT format (lon lat Z for each patch, separated by a >)
2) initialize the class changeFault
3) vary the initial fault geometry using varyDip, varyStrike or varyPosition

#### exemple:
```python
# import the classes
from epistemic_uncertainties import changeFault as chf
from epistemic_uncertainties import uncertainties as unc

# Define fault geometry parameters
faultdir = '/path/to/files'
faultfile = 'faultfile.rectangle'
nstrike = 14
utmzone = 32
dip_range = [0,1,2]
# Initialize the class changeFault
faultvar = chf(faultdir,faultfile,nstrike,utmzone,inputCoordinates='lonlatz', nbrcorners=4)
# Write new fault geometries with a dip deeper of 1 and 2 degrees
faultvar.varyDip(dip_range)
```

4) Once several fault geometries are calculated, you can calculate the related Green's functions using your favorite tool
5) Save each Green's function to a file

#### exemple:
```python
# Import you data
data = [gps]

# calculate Green's functions for each fault dip 
for d in dip_range:
    fault = rectFault('fault') 
    fault.readPatchesFromFile(faultdir+'faultfile_dip_'+str(d)+'.rectangle')   # personnal functions to import fault
    fault.buildGFs(data)  # Calculate the Green's functions
    fault.G.tofile('G_dip_'+str(d)+'.txt') # Save each GF to file
```

6) Calculate or import your a priori slip model

#### exemple:
```python
length = 26
width = 20
meanslip = -m0*10**(-7)/(3*10**10*length*width)
mprior = np.empty(len(fault.patch)*2)
mprior [0:len(fault.patch)] = 0  # strike-slip 
mprior [len(fault.patch):len(fault.patch)*2] = meanslip*100  # dip-slip, in cm!
```

Then, use the class uncertainties to calculate the uncertainties in your fault geometry:

7) Initialize the class uncertainties
8) Calculate the uncertainties with the functions calcCpDip, calcCpStrike or calcCpPosition

#### exemple:
```python
Np = 308  # number of parameters
GFs = ['G_dip_'+str(d)+'.txt' for d in dip_range]  # list of the Green's functions
sigma_dip = 5  # standard deviation on the dip parameter
uncert = uncer(faultdir, Np, multifault=False, export='./')
uncert.calcCpDip(GFs,dip_range,sigma_dip,mprior)
```

## License

Copyright © 2018-2023 Théa Ragon

This work  is free. You can  redistribute it and/or modify it under the terms of the MIT license. See the LICENSE file for more details.

This program  is free software.  It comes  without any warranty,  to the extent permitted by applicable law.
