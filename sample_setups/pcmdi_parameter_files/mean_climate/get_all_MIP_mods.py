import json
import glob
import json

pin = '/work/cmip-dyn/CMIP6/CMIP/historical/atmos/mon/rlut/'
pin = '/p/user_pub/xclim/CMIP6/CMIP/historical/atmos/mon/rlut/'

MIPS = ['CMIP5','CMIP6']
exps = ['historical','amip']

mod_dic = {}
for mip in MIPS:
 mod_dic[mip] = {}
 for exp in exps:
   mod_dic[mip][exp] = []
   ptmp = pin.replace('CMIP6',mip)
   ptmp = ptmp.replace('historical',exp)

   lst = glob.glob(ptmp +  '*.r1*.xml')  
   mods = []
   for l in lst:
#   print(l.split('.'))
    mod = l.split('.')[4]
    if mod not in mods: mods.append(mod)

   mod_dic[mip][exp] = mods

json.dump(mod_dic,open('all_mip_mods.json','w'))


  
