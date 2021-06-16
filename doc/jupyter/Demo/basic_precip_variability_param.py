mip = "cmip5"
exp = "historical"
mod = "GISS-E2-H"
var = "pr"
frq = "day"
modpath = 'demo_data_standard/CMIP5_demo_precip_var'
outdir = 'demo_output/precip_variability'
prd = [2000,2001]  # analysis period
fac = 86400  # factor to make unit of [mm/day]

# length of segment in power spectra (~10 years)
# shortened to 2 years for demo purposes
nperseg = 2 * 365
# length of overlap between segments in power spectra (~5 years)
# shortened to 1 year for demo purposes
noverlap = 1 * 365
