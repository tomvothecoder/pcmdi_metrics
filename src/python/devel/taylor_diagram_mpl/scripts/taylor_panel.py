#!/usr/bin/env python

####################################################################
# Usage example:
# python taylor_panel.py -j $(JSON_FILE) -s all -e amip -d global -o $(outputDir) -t $(tool)
# ex: python taylor_panel.py -j /work/gleckler1/processed_data/cmip5clims_metrics_package-amip/cmec_11022017/pr_2.5x2.5_regrid2_regrid2_metrics.json -s all -e amip -d global -o test2 -t uvcdat
####################################################################

import numpy as NP
import matplotlib.pyplot as PLT
import json
import sys, os
import pcmdi_metrics
from pcmdi_metrics.taylor_diagram_mpl import TaylorDiagram
import argparse
from argparse import RawTextHelpFormatter
import sys
import string
import MV2

#------------------------------------------------------------------------
P = argparse.ArgumentParser(
    description='Runs PCMDI Metrics Computations',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

P.add_argument("-j", "--json",
                      type = str,
                      dest = 'json',
                      help = "Path to json file")
P.add_argument("-v", "--variable",
                      type = str,
                      dest = 'var',
                      help = "(Case Insensitive)")
P.add_argument("-s", "--season",
                      type = str,
                      default = 'DJF',
                      help = "Season for mode of variability\n"
                             "- Available options: DJF (default), MAM, JJA, SON or all")
P.add_argument("-e", "--experiment",
                      type = str,
                      dest = 'exp',
                      default = 'historical',
                      help = "AMIP, historical or picontrol")
P.add_argument("-d", "--domain",
                      type = str,
                      dest = 'dom',
                      default = 'global',
                      help = "put options here")
P.add_argument("-o", "--plotpath",
                      type = str,
                      dest = 'outpath',
                      default = '',
                      help = "")
P.add_argument("-t", "--tool",
                      type = str,
                      dest = 'tool',
                      default = 'matplotlib',
                      help = "Available options: matplotlib, uvcdat")

args = P.parse_args(sys.argv[1:])

json_path = args.json
var = args.var
dom = args.dom
exp = args.exp
season = args.season
outpath = args.outpath 
tool = args.tool

print 'after args'
print json_path,' ',season,' ', outpath,' ', exp,' ', var , ' ', dom, ' ', tool

#------------------------------------------------------------------------
fj = open(json_path)
dd = json.loads(fj.read())
fj.close()

var = dd['Variable']['id']
#print var

mods = dd['RESULTS'].keys()
#print mods

source_ref = dd['RESULTS'][mods[0]]["defaultReference"]['source']
#print source_ref

print exp, var, json_path

#------------------------------------------------------------------------
# Seasons ---
#------------------------------------------------------------------------
if season == 'all':
  seasons = ['djf', 'mam', 'jja', 'son']
  fig_filename = var + '_' + exp + '_taylor_4panel_' + season + '_' + dom
else:
  seasons = [season]
  fig_filename = var + '_' + exp + '_taylor_1panel_' + season + '_' + dom

#------------------------------------------------------------------------
# tool ---
#------------------------------------------------------------------------
if tool == 'matplotlib':
  fig.suptitle(var.title()+', '+(exp+', '+dom).upper(), size='x-large') # Giving title for the entire canvas

  if season == 'all':
    rects = {'djf':221, 'mam':222, 'jja':223, 'son':224} # subplot location
    fig = PLT.figure(figsize=(14,11)) # optimized figure size for four subplots
  else:
    rects = {}
    rects[season] = 111 # subplot location
    fig = PLT.figure(figsize=(11,8)) # optimized figure size for one subplot

elif tool == 'uvcdat':
  import vcs
  import EzTemplate

  #canvas = vcs.init(geometry=(1200,800))
  canvas = vcs.init()
  #my_template = vcs.createtemplate()
  my_template = vcs.createtemplate(source="deftaylor")

  ## EzTemplate ---
  if season == 'all':
    M = EzTemplate.Multi(template=my_template, rows=2,columns=2, x=canvas)
  else:
    M = EzTemplate.Multi(template=my_template, rows=1,columns=1, x=canvas)

  # Legend colorbar ---
  M.legend.thickness = 0.4 # Thickness of legend color bar
  M.legend.direction = 'horizontal'

  # Border margin for entire canvas ---
  M.margins.top = .14
  M.margins.bottom = .09
  #M.margins.left = .04
  #M.margins.right = .05

#========================================================================
# Individual TD plots
#------------------------------------------------------------------------
for s, season in enumerate(seasons):
  # Reference std from obs
  stdrefs = float(dd['RESULTS'][mods[0]]["defaultReference"]['r1i1p1'][dom]['std-obs_xy'][season])

  #------------------------------------------------------------------------
  # Read in 
  #------------------------------------------------------------------------
  all_mods = []

  for mod in mods:
    cor = float(dd['RESULTS'][mod]["defaultReference"]['r1i1p1'][dom]['cor_xy'][season])
    std = float(dd['RESULTS'][mod]["defaultReference"]['r1i1p1'][dom]['std_xy'][season])
    all_mods.append([std,cor,str(mod)])

  #------------------------------------------------------------------------
  # plot 
  #------------------------------------------------------------------------
  # IF MatPlotLib
  #------------------------------------------------------------------------
  if tool == 'matplotlib':
    colors = PLT.matplotlib.cm.Set1(NP.linspace(0,1,len(all_mods)))

    dia = TaylorDiagram(stdrefs, fig=fig, rect=rects[season], label=source_ref)

    # Add samples to Taylor diagram
    for i,(stddev,corrcoef,name) in enumerate(all_mods):
      dia.add_sample(stddev, corrcoef,
                     marker='$%d$' % (i+1), ms=10, ls='',
                     #mfc='k', mec='k', # B&W
                     #mfc=colors[i], mec=colors[i], # Colors
                     label=name)

    # Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5') # 5 levels
    dia.ax.clabel(contours, inline=1, fontsize=10, fmt='%.1f')
    # Tricky: ax is the polar ax (used for plots), _ax is the
    # container (used for layout)
    dia._ax.set_title(season.upper()) # Title for the subplot

  #------------------------------------------------------------------------
  # IF UV-CDAT
  #------------------------------------------------------------------------
  elif tool == 'uvcdat':
    mods2 = list(source_ref) + mods
      
    std_list = [stdrefs]
    cor_list = [1.]

    for i,(stddev,corrcoef,name) in enumerate(all_mods):
      std_list.append(stddev)
      cor_list.append(corrcoef)

    TDdata = MV2.array(zip(std_list, cor_list)) 
    TDdata.id = var+', '+exp+', '+dom+', '+season

    if s==0:
      r=0
      c=0
    elif s==1:
      r=0
      c=1
    elif s==2:
      r=1
      c=0
    else:
      r=1
      c=1  

    #t = M.get(legend='local',row = r, column = c) # Use local colorbar
    t = M.get(legend='local',row = r, column = c) # Use local colorbar

    taylor = vcs.createtaylordiagram()
    taylor.referencevalue = stdrefs
    canvas.plot(TDdata, t, taylor, skill=taylor.defaultSkillFunction)

  #------------------------------------------------------------------------
  else:
    sys.exit(tool+' is not available option for tool (-t)')

#========================================================================
# Legend and Save
#------------------------------------------------------------------------
if not os.path.exists(outpath): os.makedirs(outpath)

#------------------------------------------------------------------------
# IF MatPlotLib
#------------------------------------------------------------------------
if tool == 'matplotlib':
  # Add a figure legend and title. For loc option, place x,y tuple inside [ ].
  # Can also use special options here:
  # http://matplotlib.sourceforge.net/users/legend_guide.html
  
  fig.legend(dia.samplePoints,
             [ p.get_label() for p in dia.samplePoints ],
             numpoints=1, prop=dict(size='small'), loc='right')
  
  PLT.savefig(outpath + '/' + fig_filename + '.png')

#------------------------------------------------------------------------
# IF UV-CDAT
#------------------------------------------------------------------------
else:
  # Logos ---
  # PCMDI
  logo2 = vcs.utils.Logo('/work/lee1043/cdat/pmp/mean_climate_maps/lib/160915_PCMDI_logo-oblong_377x300px.png')
  #logo2.x = .06
  logo2.x = .92
  logo2.y = .93
  logo2.width = logo2.source_width * .3
  logo2.height = logo2.source_height * .3
  logo2.plot(canvas)

  canvas.png(outpath + '/' + fig_filename + '.png')
