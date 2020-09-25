#!/usr/bin/python

import numpy as np
import xarray as xr
import cdms2 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from pcmdi_metrics.pcmdi.pmp_parser import PMPParser

P = PMPParser()

#P.use("--modpath")
#P.use("--modnames")
P.use("--results_dir")
#P.use("--test_data_set")
#P.use("--reference_data_path")

P.add_argument("-e", "--experiment",
               type=str,
               dest='exp',
               default=None, #'amip', #'historical',
               help="AMIP, historical or picontrol")

P.add_argument("--modpath",
               type=str,
               dest='modpath',
               default=None,
               help="path+file")

P.add_argument("--realization",
               type=str,
               dest='realization',
               default=None,
               help="path+file")

P.add_argument("--mod_name",
               type=str,
               dest='mod_name',
               default=None,
               help="path+file")

args = P.get_parameter()
exp = args.exp
rn = args.realization
mod = args.mod_name
modpath = args.modpath

missingthresh = 0.3 # threshold of missing data fraction at which a year is thrown out 

# Given one year of precip data, calculate the number of days for half of precipitation
# Ignore years with zero precip (by setting them to NaN).
# Ignore years with more than 30% missing data
def oneyear(thisyear):
    # thisyear is one year of data, (an np array) with the time variable in the leftmost dimension
    dims=thisyear.shape
    nd=dims[0]
    missingfrac = (np.sum(np.isnan(thisyear),axis=0)/nd)
    ptot=np.sum(thisyear,axis=0)
    sortandflip=-np.sort(-thisyear,axis=0)
    cum_sum=np.cumsum(sortandflip,axis=0)
    ptotnp=np.array(ptot)
    ptotnp[np.where(ptotnp == 0)]=np.nan
    pfrac = cum_sum / np.tile(ptotnp[np.newaxis,:,:],[nd,1,1])
    ndhy = np.full((dims[1],dims[2]),np.nan)
    x=np.linspace(0,nd,num=nd+1,endpoint=True)
    z=np.array([0.0])
    for ij in range(dims[1]):
        for ik in range(dims[2]):
            p=pfrac[:,ij,ik]
            y=np.concatenate([z,p])
            ndh=np.interp(0.5,y,x)
            ndhy[ij,ik]=ndh
    ndhy[np.where(missingfrac > missingthresh)] = np.nan
    return ndhy 


# Dataset - specific loop - gather an np.array of one year of data at a time, hand off to function. 

# concatenated version of files available here: ftp://ftp.climserv.ipsl.polytechnique.fr/FROGs/1DD_V0/3B42_v7.0/
# FROGS - Roca et al 2019 
# paper: https://www.earth-syst-sci-data.net/11/1017/2019/
# data: https://doi.org/10.14768/06337394-73A9-407C-9997-0E380DAC5598

diri="~/Downloads/"
diri="/p/user_pub/pmp/pmp_obs_preparation/orig/data/FROGS_precip/CMORPH_v1.0_CRT/"
diri = modpath

print('diri is ', diri)

years=[]
for year in range(1998,2011+1):
    years.append(str(year))
    
ny=len(years)
file=diri+"3B42v7.1DD."+years[0]+".nc"
file= diri+ "CMORPH_V1.0.1DD." + years[0] + ".nc"

ds = xr.open_dataset(file)
cfy = np.full((ny,ds.lon.size,ds.lat.size),np.nan)

for year in range(ny):
    print(years[year])
    file=diri+"3B42v7.1DD."+years[year]+".nc"
    file= diri+ "CMORPH_V1.0.1DD." + years[year] + ".nc"

    ds = xr.open_dataset(file)
    thisyear=ds.rain
    thisyearnp=np.array(thisyear.transpose('time','lon','lat'))
    ndhy=oneyear(thisyearnp)
    cfy[year,:,:]=ndhy    

ndm=np.nanmedian(cfy,axis=0) # ignore years with zero precip
#missingthresh = 0.3 # threshold of missing data fraction at which a location is thrown out 
missingfrac = (np.sum(np.isnan(cfy),axis=0)/ny)
ndm[np.where(missingfrac > missingthresh)] = np.nan

ndmda = xr.DataArray(ndm, dims=['lon','lat'], coords={'lon': ds.lon,'lat' : ds.lat})
ndmda.to_netcdf('nday_trmm.nc')

#plt.imshow(ndm)
plt.imshow(np.flipud(ndm.transpose()))
plt.colorbar()

filename="unevendays_trmm.pdf"
plt.savefig(filename)
plt.close()
