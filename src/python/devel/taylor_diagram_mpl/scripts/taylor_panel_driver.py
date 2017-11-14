import shlex, subprocess
import sys, os

from os import listdir
from os.path import isfile, join

for exp in ['amip', 'historical', 'picontrol']:

  mypath = '/work/gleckler1/processed_data/cmip5clims_metrics_package-'+exp+'/cmec_11022017'
  onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

  #print onlyfiles

  outpath = '/work/lee1043/cdat/pmp/clim_plots/cmip5/'+exp+'/clim/TaylorDiagram'
  if not os.path.exists(outpath): os.makedirs(outpath)

  for onlyfile in onlyfiles:
    jsonfile = os.path.join(mypath, onlyfile)
    
    command_line = 'python taylor_panel.py -j '+jsonfile+' -s all -e '+exp+' -d global -o '+outpath
    args = shlex.split(command_line)

    p = subprocess.Popen(args)
