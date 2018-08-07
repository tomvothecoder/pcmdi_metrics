#!/usr/bin/env python

#########################################################
# SAMPLE COMMAND LINE EXECUTION USING ARGUMENTS BELOW
#########################################################
# python cloud_metrics_zelinka2013.py 
# -mp /work/cmip5/piControl/atm/mo/clisccp/cmip5.MODS.piControl.r1i1p1.mo.atm.cfMon.clisccp.ver-1.latestX.xml
# -op /work/clim_obs/orig/data/ISCCP/clisccp_198307-200806.nc
# -kp /work/zhang24/pmp_cloud_K13/data/cloud_kernels2.nc
# --mns MIROC5 MIROC-ESM-CHEM
# --var clisccp
# --varobs clisccp (varobs needed only when varname is different to model in obs) 
# --outpd /work/zhang24/pmp_cloud_K13/test
# --outpj /work/zhang24/pmp_cloud_K13/test 
# --outnj output.json 
#########################################################

import logging
LOG_LEVEL = logging.INFO
logging.basicConfig(level=LOG_LEVEL)

import cdms2 as cdms
import copy
import sys
import os
import string
import json
import pcmdi_metrics
from pcmdi_metrics.pcmdi.pmp_parser import PMPParser
import collections
from collections import defaultdict

import cdutil
import MV2 as MV
import numpy as np
import pylab as pl



#debug = True
debug = False

def tree(): return defaultdict(tree)

#########################################################

P = PMPParser() # Includes all default options

P.add_argument("--mp", "--modpath",
               type=str,
               dest='modpath',
               required=True,
               help="Explicit path to model monthly clisccp time series")
P.add_argument("--op", "--obspath",
               type=str,
               dest='obspath',
               default='',
               help="Explicit path to obs monthly clisccp time series")
P.add_argument("--kp", "--kernelpath",
               type=str,
               dest='kernelpath',
               default='',
               help="Explicit path to the Zelinka et al 2012 kernels")
P.add_argument('--mns', '--modnames',
               type=str,
               nargs='+',
               dest='modnames',
               required=True,
               help='Models to apply')
P.add_argument("--var", "--variable",
               type=str,
               dest='variable',
               default='clisccp',
               help="Variable: 'clisccp'")
P.add_argument("--varobs", "--variableobs",
               type=str,
               dest='variableobs',
               default='',
               help="Variable name in observation (default: same as var)")
P.add_argument("--outpj", "--outpathjsons",
               type=str,
               dest='outpathjsons',
               default='.',
               help="Output path for jsons")
P.add_argument("--outnj", "--outnamejson",
               type=str,
               dest='jsonname',
               default='cloud_K13.json',
               help="Output path for jsons")
P.add_argument("--outpd", "--outpathdata",
               type=str,
               dest='outpathdata',
               default='.',
               help="Output path for data")
P.add_argument("-e", "--experiment",
               type=str,
               dest='experiment',
               default='historical',
               help="AMIP, historical or picontrol")
P.add_argument("-c", "--MIP",
               type=str,
               dest='mip',
               default='CMIP5',
               help="put options here")
P.add_argument("-p", "--parameters",
               type=str,
               dest='parameters',
               default='',
               help="")

args = P.parse_args(sys.argv[1:])

modpath = args.modpath
obspath = args.obspath
kernelpath = args.kernelpath
mods = args.modnames
var = args.variable
varobs = args.variableobs
if varobs == '': varobs = var
outpathjsons = args.outpathjsons
outfilejson = args.jsonname
outpathdata = args.outpathdata
exp = args.experiment

###########################################################################
# HELPFUL FUNCTIONS FOLLOW 
###########################################################################

###########################################################################
def add_cyclic(data):
    # Add Cyclic point around 360 degrees longitude:

    # This function assumes that your longitudes range from 0 to 360, not -180 to 180
    lons=data.getLongitude()[:]
    dx=np.gradient(lons)[-1]
    data2 = data(longitude=(0, dx+np.max(lons)), squeeze=True)    
    return data2

###########################################################################
def reshape_generic(orig_data,data_to_match):

    # this function resizes and tiles orig_data the same shape as data_to_match
    # orig_data must have fewer dimensions than data_to_match

    A=orig_data.shape
    B=data_to_match.shape
    ndim_new = data_to_match.ndim

    # find index where B disagrees with A
    #disagree=np.setdiff1d(B,A)
    agree=np.in1d(B,A)
    j=[]
    for i in range(len(B)):
        ndim_orig = orig_data.ndim
        if agree[i]==False:
            j=np.append(j,i)
            new = np.expand_dims(orig_data,axis=ndim_orig)
            NEW = np.tile(new,[B[i]])
            new_mask = np.expand_dims(orig_data.mask,axis=ndim_orig)
            MASK = np.tile(new_mask,[B[i]])
            orig_data = np.ma.array(NEW,mask=MASK)

    # need to move axes around
    for i in range(len(B)):
        C=orig_data.shape
        if C[i]!=B[i]:
            orig_data = np.moveaxis(orig_data, i, B.index(C[i])) # (a, source, destination)

    return orig_data

###########################################################################
def nanarray(vector):

    # this generates a masked array with the size given by vector
    # example: vector = (90,144,28)

    # similar to this=NaN*ones(x,y,z) in matlab

    this=MV.zeros(vector)
    this=MV.masked_where(this==0,this)

    return this

###########################################################################
def map_SWkern_to_lon(Ksw,albcsmap):

    from scipy.interpolate import interp1d
    ## Map each location's clear-sky surface albedo to the correct albedo bin
    # Ksw is size 12,7,7,lats,3
    # albcsmap is size A,lats,lons
    albcs=np.arange(0.0,1.5,0.5) 
    A=albcsmap.shape[0]
    TT=Ksw.shape[1]
    PP=Ksw.shape[2]
    lenlat=Ksw.shape[3]
    lenlon=albcsmap.shape[2]
    SWkernel_map=nanarray((A,TT,PP,lenlat,lenlon))
    for M in range(A):
        MM=M
        while MM>11:
            MM=MM-12
        for LA in range(lenlat):
            alon=albcsmap[M,LA,:] 
            # interp1d can't handle mask but it can deal with NaN (?)
            try:
                alon2=MV.where(alon.mask,np.nan,alon)   
            except:
                alon2=alon
            if np.ma.count(alon2)>1: # at least 1 unmasked value
                if len(pl.find(Ksw[MM,:,:,LA,:]>0))==0:
                    SWkernel_map[M,:,:,LA,:] = 0
                else:
                    f = interp1d(albcs,Ksw[MM,:,:,LA,:],axis=2)
                    ynew = f(alon2.data)
                    ynew=MV.masked_where(alon2.mask,ynew)
                    SWkernel_map[M,:,:,LA,:] = ynew
            else:
                continue

    return SWkernel_map

##########################################################
##########################################################

if var != 'clisccp' :
    sys.exit('Variable '+var+' is not correct')

# SETUP WHERE TO OUTPUT RESULTING  (netcdf)
try:
    jout = outpathjsons
    os.mkdir(jout)
except BaseException:
    pass

models = copy.copy(args.modnames)
print models
#if obspath != '':
#models.insert(0,'obs')

# DICTIONARY TO SAVE RESULT
cloud_metrics_dic = tree() ## Set tree structure dictionary

#=================================================
# Loop for Observation and Models 
#-------------------------------------------------
for mod in models:
    print ' ----- ', mod,' ---------------------'
  
    if mod == 'obs':
        file_path = obspath
        varname = varobs
        mods_key = 'OBSERVATION'
    else:
        file_path = modpath.replace('MODS', mod)
        file_path_rsdscs = file_path.replace('clisccp','rsdscs')
        file_path_rsdscs = file_path_rsdscs.replace('cfMon','Amon')
        file_path_rsuscs = file_path.replace('clisccp','rsuscs')
        file_path_rsuscs = file_path_rsuscs.replace('cfMon','Amon')
        varname = var
        mods_key = 'MODELS'
  
    try:
  
        #f = cdms.open(file_path)
        print 'test: ',file_path
        cloud_metrics_dic[mods_key][mod]['input_data'] = file_path
    
        if debug: print file_path 

###########################################################################
# MAIN ROUTINE FOLLOWS
###########################################################################

        # Get Kernel Grid (for purpose of regridding):
        #=====================
        print 'entered test2!'
        kdir='/work/zhang24/pmp_cloud_K13/data/sw_a_gfdl_std.nc'
        f=cdms.open(kdir,'r')
        kern_grid = f('sw_a').getGrid()
        f.close()

        # Load in the Zelinka et al 2012 kernels:
        f=cdms.open(kernelpath)
        LWkernel0=f('LWkernel')
        SWkernel0=f('SWkernel')
        f.close()
        print 'test2: ',kern_grid.shape,SWkernel0.shape

        # Take only the portion of the kernel histogram where there are obs (ignore first tau bin)
        SWkernel = SWkernel0[:,1:,:]
        LWkernel = LWkernel0[:,1:,:]
        del(LWkernel0,SWkernel0)

        albcs=np.arange(0.0,1.5,0.5) # the clear-sky albedos over which the kernel is computed

        ######################################################
        ############# Load in ISCCP observations #############
        ######################################################
        Itslc=("1984-01-01","2007-12-31") # Full years of overlap between model and ISCCP
        f=cdms.open(obspath,'r')
        obs_clisccp0=f(varname,longitude=(0,359.9),time=Itslc,squeeze=1)
        f.close()
        print 'test3: ',obs_clisccp0.shape

        # Compute Climatological Annual Cycle:
        cdutil.setTimeBoundsMonthly(obs_clisccp0)
        obs_clisccp = cdutil.ANNUALCYCLE.climatology(obs_clisccp0) #(12,...)
        del(obs_clisccp0)

        # Regrid and flip the CTP dimension to go SFC to TOA:
        obs_clisccp = add_cyclic(obs_clisccp)
        obs_clisccp_grd = obs_clisccp[:,:,-1::-1,:].regrid(kern_grid,regridTool="esmf",regridMethod = "linear")


        agg_mod_clisccp_bias = nanarray((12,2,3,60,144)) # (month, tau_bins, CTP_bins, lat, lon)
        agg_mod_SW_bias = nanarray((12,2,3,60,144)) # (month, tau_bins, CTP_bins, lat, lon)
        agg_mod_LW_bias = nanarray((12,2,3,60,144)) # (month, tau_bins, CTP_bins, lat, lon)

        agg_obs_clisccp_bias=nanarray((12,2,3,60,144))
        agg_obs_SW_bias=nanarray((12,2,3,60,144))
        agg_obs_LW_bias=nanarray((12,2,3,60,144))

        ######################################################
        ############# Load in CLISCCP from model #############
        ######################################################
        # Grab a random AMIP simulation:
        print 'enter test4! ',file_path,varname
        f=cdms.open(file_path,'r')
        clisccp0=f(varname,time=Itslc,squeeze=1)
        f.close()
        print 'test4: ',clisccp0.shape

        # Compute Climatological Annual Cycle:
        clisccp = cdutil.ANNUALCYCLE.climatology(clisccp0) #(12,...)
        del(clisccp0)

        # Remove the thinnest optical depth bin so as to compare properly with obs:
        clisccp=clisccp[:,1:,:,:]

        # Make sure cloud fractions are in percent  
        sumclisccp=MV.sum(MV.sum(clisccp,2),1)
        if np.max(sumclisccp) <= 1.:
            clisccp = clisccp*100.          

        ######################################################
        ########## Compute clear-sky surface albedo ##########
        ######################################################
        print 'enter test5! ',file_path_rsdscs
        f=cdms.open(file_path_rsdscs,'r')
        rsdscs0 = f('rsdscs',squeeze=1) # Clearsky downwelling solar flux at surface
        f.close()
        print 'middle test5!'
        f=cdms.open(file_path_rsuscs,'r')
        rsuscs0 = f('rsuscs',squeeze=1) # Clearsky upwelling solar flux at surface
        f.close()
        print 'test5: ',rsdscs0.shape

        # Compute Climatological Annual Cycle:
        rsdscs = cdutil.ANNUALCYCLE.climatology(rsdscs0) #(12,...)
        rsuscs = cdutil.ANNUALCYCLE.climatology(rsuscs0) #(12,...)

        albcs = rsuscs/rsdscs
        albcs=MV.where(albcs>1.,1,albcs) # where(condition, x, y) is x where condition is true, y otherwise
        albcs=MV.where(albcs<0.,0,albcs)

 
        # Regrid everything to the kernel grid:
        albcs = add_cyclic(albcs)
        albcs_grd = albcs.regrid(kern_grid,regridTool="esmf",regridMethod = "linear")
        clisccp = add_cyclic(clisccp)
        clisccp_grd = clisccp.regrid(kern_grid,regridTool="esmf",regridMethod = "linear")

        ## Use average control albcs to map SW kernel to appropriate longitudes
        SWkernel_map = map_SWkern_to_lon(SWkernel,albcs_grd)

        # LW kernel does not depend on albcs, just repeat the final dimension over longitudes:
        A=SWkernel_map.shape[0]
        LWkernel_map0=np.tile(np.tile(LWkernel[:,:,:,:,0],(1,1,1,1,1)),(144,1,1,1,1))(order=[1,2,3,4,0])
        LWkernel_map=nanarray(SWkernel_map.shape)
        for a in range(A):
            aa=a
            while aa>11:
                aa=aa-12
            LWkernel_map[a,:] = LWkernel_map0[aa,:]

        ## Compute Cloud Fraction Histogram Anomalies w.r.t. observations
        clisccp_bias = clisccp_grd - obs_clisccp_grd

        ## Multiply Anomalies by Kernels
        SW0 = SWkernel_map*clisccp_bias
        LW = LWkernel_map*clisccp_bias

        ## Set the SW cloud feedbacks to zero in the polar night
        # The sun is down if every bin of the SW kernel is zero:
        sundown=MV.sum(MV.sum(SWkernel_map,axis=2),axis=1)  #MO,90,144
        repsundown=np.tile(np.tile(sundown,(1,1,1,1,1)),(7,6,1,1,1))(order=[2,1,0,3,4])
        SW1 = MV.where(repsundown==0, 0, SW0) # where(condition, x, y) is x where condition is true, y otherwise
        SW = MV.where(repsundown.mask, 0, SW1) # where(condition, x, y) is x where condition is true, y otherwise

        # SW and LW contain the SW and LW radiation anomalies contributed from cloud anomalies in each bin of the histogram
        LW.setAxisList(clisccp_bias.getAxisList())
        SW.setAxisList(clisccp_bias.getAxisList())
        LWkernel_map.setAxisList(clisccp_bias.getAxisList())
        SWkernel_map.setAxisList(clisccp_bias.getAxisList())


        ########################################################
        ######### Compute Klein et al (2013) metrics ###########
        ########################################################
        eq60 = cdutil.region.domain(latitude=(-60.,60.)) # equatorward of 60

        ## E_TCA (TOTAL CLOUD AMOUNT METRIC)
        # take only clouds with tau>1.3 between 60S-60N
        obs_clisccp_eq60 = eq60.select(obs_clisccp_grd[:,1:,:])
        mod_clisccp_eq60 = eq60.select(clisccp_grd[:,1:,:])

        # sum over CTP and TAU:
        mod_cltisccp_eq60 = MV.sum(MV.sum(mod_clisccp_eq60,axis=1),axis=1) # (time, lat, lon)
        obs_cltisccp_eq60 = MV.sum(MV.sum(obs_clisccp_eq60,axis=1),axis=1) # (time, lat, lon)
        obs_cltisccp_eq60 = MV.sum(MV.sum(obs_clisccp_eq60,axis=1),axis=1) # (time, lat, lon)

        ########################################################
        # E_TCA for Model minus ISCCP:
        ########################################################
        # 1) Denominator (Eq. 3 in Klein et al. (2013))
        avg = cdutil.averager(MV.average(obs_cltisccp_eq60,axis=0), axis='xy', weights='weighted') # (scalar)
        rep_avg = reshape_generic(avg,obs_cltisccp_eq60) # (time, lat, lon)
        anom = obs_cltisccp_eq60 - rep_avg # anomaly of obs from its spatio-temporal mean
        E_TCA_denom = np.ma.sqrt(cdutil.averager(MV.average(anom**2,axis=0), axis='xy', weights='weighted')) # (scalar)
        # 2) Numerator
        anom = mod_cltisccp_eq60 - obs_cltisccp_eq60  # (time, lat, lon)
        E_TCA_numer = np.ma.sqrt(cdutil.averager(MV.average(anom**2,axis=0), axis='xy', weights='weighted')) # (scalar)
        E_TCA_mod = E_TCA_numer/E_TCA_denom
     
        # E_TCA for MODIS minus ISCCP (where they overlap):

        ########################################################
        # CLOUD PROPERTY METRICS
        ########################################################
        # take only clouds with tau>3.6 between 60S-60N
        clisccp_bias_eq60 = eq60.select(clisccp_bias[:,2:,:])
        obs_clisccp_eq60 = eq60.select(obs_clisccp_grd[:,2:,:])
        mod_clisccp_eq60 = eq60.select(clisccp_grd[:,2:,:])
        LWkernel_eq60 = eq60.select(LWkernel_map[:,2:,:])
        SWkernel_eq60 = eq60.select(SWkernel_map[:,2:,:])

        # Compute anomaly of obs histogram from its spatio-temporal mean
        avg_obs_clisccp_eq60 = cdutil.averager(obs_clisccp_eq60, axis='xy', weights='weighted') # (time,TAU,CTP)
        rep_avg_obs_clisccp_eq60 = reshape_generic(avg_obs_clisccp_eq60,obs_clisccp_eq60) # (time, TAU, CTP, lat, lon)
        anom_obs_clisccp_eq60 = obs_clisccp_eq60 - rep_avg_obs_clisccp_eq60 # anomaly of obs from its spatio-temporal mean

        ## Compute radiative impacts of cloud fraction anomalies
        mod_SW_bias = eq60.select(SW[:,2:,:])
        obs_SW_bias = anom_obs_clisccp_eq60*SWkernel_eq60
        mod_LW_bias = eq60.select(LW[:,2:,:])
        obs_LW_bias = anom_obs_clisccp_eq60*LWkernel_eq60

        ## Aggregate high, mid, and low clouds over medium and thick ISCCP ranges
        Psec_name = ['LO','MID','HI']
        Psections=[slice(0,2),slice(2,4),slice(4,7)]
        Psec_dic=dict(zip(Psec_name,Psections))
        Tsec_name = ['MED','THICK']
        Tsections=[slice(0,2),slice(2,4)]
        Tsec_dic=dict(zip(Tsec_name,Tsections))

        tt=-1
        for Tsec in Tsec_name:
            tt+=1
            TT=Tsec_dic[Tsec]
            pp=-1
            for Psec in Psec_name:
                pp+=1
                PP=Psec_dic[Psec]
                agg_obs_SW_bias[:,tt,pp,:] = MV.sum(MV.sum(obs_SW_bias[:,TT,PP,:],axis=1),axis=1)
                agg_mod_SW_bias[:,tt,pp,:] = MV.sum(MV.sum(mod_SW_bias[:,TT,PP,:],axis=1),axis=1)
                agg_obs_LW_bias[:,tt,pp,:] = MV.sum(MV.sum(obs_LW_bias[:,TT,PP,:],axis=1),axis=1)
                agg_mod_LW_bias[:,tt,pp,:] = MV.sum(MV.sum(mod_LW_bias[:,TT,PP,:],axis=1),axis=1)
                agg_obs_clisccp_bias[:,tt,pp,:] = MV.sum(MV.sum(anom_obs_clisccp_eq60[:,TT,PP,:],axis=1),axis=1)
                agg_mod_clisccp_bias[:,tt,pp,:] = MV.sum(MV.sum(clisccp_bias_eq60[:,TT,PP,:],axis=1),axis=1)

        ## Compute E_ctp-tau -- Cloud properties error
        ctot = MV.sum(MV.sum(agg_mod_clisccp_bias**2,axis=1),axis=1)/6;
        ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
        E_ctpt_numer = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

        ctot = MV.sum(MV.sum(agg_obs_clisccp_bias**2,axis=1),axis=1)/6;
        ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
        E_ctpt_denom = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

        E_ctpt_mod = E_ctpt_numer/E_ctpt_denom

        ## Compute E_LW -- LW-relevant cloud properties error
        ctot = MV.sum(MV.sum(agg_mod_LW_bias**2,axis=1),axis=1)/6;
        ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
        E_LW_numer = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

        ctot = MV.sum(MV.sum(agg_obs_LW_bias**2,axis=1),axis=1)/6;
        ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
        E_LW_denom = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

        E_LW_mod = E_LW_numer/E_LW_denom

        ## Compute E_SW -- SW-relevant cloud properties error
        ctot = MV.sum(MV.sum(agg_mod_SW_bias**2,axis=1),axis=1)/6;
        ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
        E_SW_numer = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

        ctot = MV.sum(MV.sum(agg_obs_SW_bias**2,axis=1),axis=1)/6;
        ctot.setAxisList(clisccp_bias_eq60[:,0,0,:].getAxisList())
        E_SW_denom = np.ma.sqrt(cdutil.averager(ctot, axis='xy', weights='weighted')) # (time)

        E_SW_mod = E_SW_numer/E_SW_denom

      
      
        # Record cloud metrics from above calculation ---
        cloud_metrics_dic[mods_key][mod]['E_TCA'] = float(E_TCA_mod)
        cloud_metrics_dic[mods_key][mod]['E_CP'] = E_ctpt_mod.tolist()
        cloud_metrics_dic[mods_key][mod]['E_CP_axis'] = str(E_ctpt_mod.getAxis(0))
        cloud_metrics_dic[mods_key][mod]['E_LWCP'] = E_LW_mod.tolist()
        cloud_metrics_dic[mods_key][mod]['E_LWCP_axis'] = str(E_LW_mod.getAxis(0))
        cloud_metrics_dic[mods_key][mod]['E_SWCP'] = E_SW_mod.tolist()
        cloud_metrics_dic[mods_key][mod]['E_SWCP_axis'] = str(E_SW_mod.getAxis(0))
        #print 'results: ',E_TCA_mod,E_ctpt_mod,E_LW_mod,E_SW_mod
        #stop

        #f.close()
    
    except:
        print 'failed for ', mod

#=================================================
#  OUTPUT METRICS TO JSON FILE
#-------------------------------------------------
OUT = pcmdi_metrics.io.base.Base(os.path.abspath(jout), outfilejson)

disclaimer = open(
    os.path.join(
        sys.prefix,
        "share",
        "pmp",
        "disclaimer.txt")).read()

metrics_dictionary = collections.OrderedDict()
metrics_dictionary["DISCLAIMER"] = disclaimer
metrics_dictionary["REFERENCE"] = "The statistics in this file are based on Klein, et al. JGR-Atmos (2013), vol.118, 1329-1342. doi:10.1002/jgrd.50141"
metrics_dictionary["RESULTS"] = cloud_metrics_dic  # collections.OrderedDict()
print metrics_dictionary
OUT.var = var
OUT.write(
    metrics_dictionary,
    json_structure=["model", "statistic"],
    indent=3,
    separators=(
        ',',
        ': '),
    sort_keys=True)

sys.exit('done')
