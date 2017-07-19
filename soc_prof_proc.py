#!/bin/python

import numpy as np
import pandas as pd
import soc_analysis_lib as soca
import isamcalc_lib as isam
import socplot_lib as socplt
import auxiliary_lib as au

# Open and read in data
filename='N_circum_qaed.xlsx'
data = pd.read_excel(filename,index_col='Profile ID', skiprows=[1])
all_profid = data.index.unique()
lons = soca.getvarxls(data,'Long',all_profid,0)
lats = soca.getvarxls(data,'Lat N',all_profid,0)

# Sort the lons to avoid the bug from matplot
lon_s = lons[lons.argsort()]
lat_s = lats[lons.argsort()]
# Plot profile locations on the map
tit = 'Location of all SOC samples'
path = 'Location_map.png'
status = socplt.plot_profmap(-180, -60, 180, 80, lon_s, lat_s, tit, path)

# Store coordinates into ascii file, needed by GRASS
fname = './allsites.asc'
pp = np.column_stack((lons,lats))
np.savetxt(fname, pp, fmt='%8.8f', delimiter=' ', newline='\n', header='', footer='')
## Use Grass to get the permafrost status of each sites
status = au.exec_shell('bash', './get_pf_status.sh', [])

# Get the permafrost status of all sites
pfst_table = np.loadtxt(fname='N_circum_sites_pf', dtype='float', delimiter=' ')
pfst = pfst_table[:,3].astype('int')
# Only select sites with continuous and discontinuous permafrost
# 1, 2, 5, 6, 9, 10, 13, 14, 17, 18.
lgc = np.ones((len(pfst)), dtype=bool)
for i in range(0, len(pfst)):
    lgc[i] = pfst[i] in [1, 2, 5, 6, 9, 10, 13, 14, 17, 18]
sel_profid = pfst_table[:,2].astype('int')[lgc]
data_sel = soca.subset_by_ids(data, sel_profid)
# Replace space with underscore
data_sel.columns = au.rep_sp_by_sym(data_sel, sym='_')

# Extract vegetation information
sel_profid = data_sel.index.unique()
pfst_veg = np.loadtxt(fname='N_circum_sites_veg', dtype='float', delimiter=',')
site_veg = pfst_veg[:,3].astype('int')
# Transfer into pd dataframe
prof_veg = pd.DataFrame(np.stack((all_profid, site_veg), axis=1), index=all_profid, columns=['all_profid', 'veg_id'])
sel_veg = soca.subset_by_ids(prof_veg, sel_profid)

# Extract data fields for mass preserving interpolation
# We only want Lat, Lon, basal depth, layer depth, soil order, 
# suborder, Horizon type, bulk density and C density 
data_out = data_sel[['Veg_Class', 'Lat_N', 'Long', 'Basal_depth', 'Layer_depth', 'Soil_Order', \
                                'Suborder', 'Horizon_type', 'bulk_density', 'C_Density']]
# Write out CSV file, needed by R script
data_out.to_csv('./sel_sites_carbon.csv')
# Then call R script to run the mass-preseving spline interpolation
status = au.exec_shell('Rscript', 'SOCinterp.r', [])

# Read in the interpolated SOC profile
filename = 'SOCprofile.csv'
soc_interped = pd.read_csv(filename,encoding='iso-8859-1',index_col=0)
pid = soc_interped.index
# Draw SOC for each profile
soc_interped = soc_interped * 1000.   # gcm-3 to kgm-3
soc_interped.index = sel_profid

# Plot individual profile
for i in range(len(sel_profid)):
    X1 = soc_interped.loc[pid[i],:]   # SOC profile kgCm-3
    Y = np.arange(1,201)     # 1cm to 200cm
    if (data_sel.Profile_Name[sel_profid[i]].__class__.__name__ == 'Series'):
        tit = data_sel.Profile_Name[sel_profid[i]].unique().astype('string')[0]
    else:
        tit = data_sel.Profile_Name[sel_profid[i]].encode('ascii','ignore')
    path = './Figs/'+str(sel_profid[i])+'_'+tit+'.png'
    status = socplt.plot_soilprof(X1, Y, tit, path)

# Separate different sites based on ecoregion type and landcover type
# 1. Separate into different soil order + suborder
# Orthel
data_orthel, orthel_profid = soca.subset_by_cols(data_sel, 'Suborder', 'Orthel', True)
soc_orthel = soca.subset_by_ids(soc_interped, orthel_profid)
# Turbel
data_turbel, turbel_profid = soca.subset_by_cols(data_sel, 'Suborder', 'Turbel', True)
soc_turbel = soca.subset_by_ids(soc_interped, turbel_profid)
# Histel
data_histel, histel_profid = soca.subset_by_cols(data_sel, 'Suborder', 'Histel', True)
soc_histel = soca.subset_by_ids(soc_interped, histel_profid)
# Others
data_others, others_profid = soca.subset_by_cols(data_sel, 'Soil_Order', 'Gelisol', False)
soc_others = soca.subset_by_ids(soc_interped, others_profid)

# Derive mean and std.
m_orthel = soc_orthel.mean(axis=0)
s_orthel = soc_orthel.std(axis=0)
m_turbel = soc_turbel.mean(axis=0)
s_turbel = soc_turbel.std(axis=0)
m_histel = soc_histel.mean(axis=0)
s_histel = soc_histel.std(axis=0)
m_others = soc_others.mean(axis=0)
s_others = soc_others.std(axis=0)

# Make figure for each different soil order / suborder
tit = 'Orthel'
path = './Figs/'+tit+'.png'
status = socplt.plot_prof_with_errbar(m_orthel[0:199], np.arange(1,200), s_orthel[0:199], tit, path)
tit = 'Turbel'
path = './Figs/'+tit+'.png'
status = socplt.plot_prof_with_errbar(m_turbel[0:199], np.arange(1,200), s_turbel[0:199], tit, path)
tit = 'Histel'
path = './Figs/'+tit+'.png'
status = socplt.plot_prof_with_errbar(m_histel[0:199], np.arange(1,200), s_histel[0:199], tit, path)
tit = 'Others'
path = './Figs/'+tit+'.png'
status = socplt.plot_prof_with_errbar(m_others[0:199], np.arange(1,200), s_others[0:199], tit, path)

## 2. Separate into different vegetation cover under each soil order + suborder
## For orthel
## Herbacious
# lgc = np.ones((len(soc_orthel.index)), dtype=bool)
# for i in range(len(soc_orthel.index)):
#     lgc[i] = sel_profid[i] in [20]
# soc_others = soc_interped[lgc]

# 3. Randomly pick sites from each category and then use these samples to calibrate model
# Get the required information for driving the model
lons = soca.getvarxls(data_sel,'Long',sel_profid,0)
lats = soca.getvarxls(data_sel,'Lat_N',sel_profid,0)
data_sel_meta=np.column_stack((sel_profid,lons,lats))
fname='site_latlon.asc'
np.savetxt(fname, data_sel_meta, fmt='%8.8f', delimiter=' ', newline='\n', header='', footer='')

# 4. Evaluate model's performance by running all observed soil profiles and compare under each different categories
# Read in the model output
socm = pd.read_table('isam_soc.dat', header=None, delimiter=r"\s+")
socm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
socm = socm.set_index('ID')
mod_profid = socm.index
# Unit conversion from kgCm-2 to kgCm-3
z, dz, zsoih = isam.get_isam_soildp(10)
socm = socm / dz

# Get the observed profiles by using IDs from model
soco_m = np.zeros(shape=(len(mod_profid),200) , dtype=float)
for i in range(len(mod_profid)):
    soco_m[i,:] = soc_interped[sel_profid == mod_profid[i]].astype(float)
# Transfer into pandas dataframe
soco = pd.DataFrame(soco_m, index=mod_profid, columns=soc_interped.columns)
# Generate figures
for i in range(len(mod_profid)):
    Xobs = soco.iloc[i,:]   # SOC profile kgCm-3
    Yobs = np.arange(1,201)     # 1cm to 200cm
    Xmod = socm.iloc[i,:]   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    if (data_sel.Profile_Name[mod_profid[i]].__class__.__name__ == 'Series'):
        tit = data_sel.Profile_Name[mod_profid[i]].unique().astype('string')[0]
    else:
        tit = data_sel.Profile_Name[mod_profid[i]].encode('ascii','ignore')
    path = './Figs_obsvsmod/'+str(mod_profid[i])+'_'+tit+'.png'
    status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path)

# Prepare the figure showing the model results plus the observation

# 5. Calculate the mean and std then make the plot.

