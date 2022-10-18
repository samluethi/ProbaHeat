#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 17 2022
@author: Samuel Luethi - samuel.luethi@usys.ethz.ch

"""
import numpy as np
import pandas as pd
import pickle
import os

from pathlib import Path

''' load & prepare data '''
# load obs datar
DP = Path('/Users/sam/OneDrive - ETH Zurich/WCR/Projects/Heat/Data/')
TP = DP.joinpath('ESM/Global/run_220322/')
all_mort = pd.ExcelFile(DP.joinpath('MCC/clean_data/all_MORT_run_220328.xlsx'))
all_rr = pd.ExcelFile(DP.joinpath('MCC/clean_data/all_RR_run_220328.xlsx'))
all_tmm = pd.read_csv(DP.joinpath('MCC/clean_data/all_TMM_220328.csv'))
cities = all_mort.sheet_names # reduced
city_meta = pd.read_excel(DP.joinpath('MCC/MCC_Metadata_20211007.xlsx'),
                          sheet_name="cities info")
city_meta = city_meta[city_meta['citiescode'].isin(cities)]
        
cityIndex = dict(zip(cities, range(len(cities))))
city_meta['city_sort'] = city_meta['citiescode'].map(cityIndex)
city_meta.sort_values(['city_sort'], ascending = True, inplace = True)
city_meta.drop('city_sort', 1, inplace = True)

all_tmm['city_sort'] = all_tmm['city'].map(cityIndex)
all_tmm.sort_values(['city_sort'], ascending = True, inplace = True)
all_tmm.drop('city_sort', 1, inplace = True)


warming_levels = [0.7, 1.2, 1.5, 2.0, 2.5]
LE_models = ['CESM12-LE', 'CESM1-CAM5', 'CanESM2', 'GFDL-ESM2M', 'CSIRO-Mk3-6-0']

'''------------------------------------------------------------------------'''
''' calc heat mortality from bias corrected temperature data '''

# calculate heat related mortality from bias corrected large ensemble data
impacts = dict()
return_period = dict()

for model in LE_models:
    print(model)
    file_path = TP.joinpath('city_temp_data_BC_'+model+'.h5')
    imp_m = dict()
    
    for cc, city in enumerate(cities):
        if np.mod(cc, 25)==0: print(cc)
        RR = pd.read_excel(all_rr, sheet_name=city)
        mort = pd.read_excel(all_mort, sheet_name=city)
        imp_c = dict()
        
        for w, wl in enumerate(warming_levels):
            T = pd.read_hdf(file_path, key=city+'/period_'+str(w))
            if model=='EC-EARTH':
                # remove leap days for EC-EARTH
                date_series = pd.to_datetime(T.columns)
                idx_leap_day = np.where((date_series.is_leap_year==True) & \
                                        (date_series.day_of_year==60))[0]
                T = T.drop(T.columns[idx_leap_day], axis=1)
            imp, rp = calc_heat_mortality(T, mort, RR, all_tmm.TMM.values[all_tmm.city==city][0])
            imp_c[str(wl)] = imp
        imp_m[city] = imp_c
    impacts[model] = imp_m
    return_period[model] = rp

# calculate changes
RP_changes = dict()
Pct_mort = dict()
Imp_changes_rel = dict()
Imp_changes_pct = dict()

for model in LE_models:
    print(model)
    rp = return_period[model]
    rpc_m = dict()
    pm_m = dict()
    icr_m = dict()
    icp_m = dict()
    for w, wl in enumerate(warming_levels):
        rpc_w = np.zeros(len(cities))
        pm_w = np.zeros(len(cities))
        icr_w = np.zeros(len(cities))
        icp_w = np.zeros(len(cities))
        for cc, city in enumerate(cities):
            mort = pd.read_excel(all_mort, sheet_name=city)
            imp_base = impacts[model][city][str(warming_levels[0])][0,:]
            imp_per = impacts[model][city][str(wl)][0,:]
            rpc_w[cc] = change_in_RP(100, rp, imp_base, rp, imp_per)
            pm_w[cc] = pct_heat_mort(100, imp_per, rp, mort)
            icr_w[cc] = change_in_mort(100, imp_base, imp_per, rp)
            icp_w[cc] = change_in_mort_pct(100, imp_base, imp_per, rp, mort)
        rpc_m[wl] = rpc_w
        pm_m[wl] = pm_w
        icr_m[wl] = icr_w
        icp_m[wl] = icp_w
    RP_changes[model] = rpc_m
    Pct_mort[model] = pm_m
    Imp_changes_rel[model] = icr_m
    Imp_changes_pct[model] = icp_m


# observational data
imp_obs = dict()
rp_obs = dict()

for cc, city in enumerate(cities):
    if np.mod(cc, 25)==0: print(cc)
    RR = pd.read_excel(all_rr, sheet_name=city)
    mort = pd.read_excel(all_mort, sheet_name=city)
    imp_o, rp_o = calc_heat_mortality_obs(mort, RR,
                                all_tmm.TMM.values[all_tmm.city==city][0])
    imp_obs[city] = imp_o
    rp_obs[city] = rp_o
    
# add median
def add_median_imp(impacts, return_period):
    # data structure
    models = list(impacts.keys())
    locs = list(impacts[models[0]].keys())
    warming_levels = list(impacts[models[0]][locs[0]].keys())
    # take return periods of ensemble with most model runds
    ens = {key: len(value) for key, value in return_period.items()}
    max_rp_model = max(ens, key=ens.get)
    rp_main = return_period[max_rp_model]
    imp_loc = dict()
    # calc median on return period
    for loc in locs:
        imp_w = dict()
        for wl in warming_levels:
            med = np.empty([len(rp_main), len(models)])
            med_low = np.empty([len(rp_main), len(models)])
            med_high = np.empty([len(rp_main), len(models)])
            for m, model in enumerate(models):
                rp = return_period[model]
                interp_imp = np.interp(rp_main, rp, impacts[model][loc][wl][0])
                interp_imp[rp_main>np.max(rp)] = np.nan
                med[:,m] = interp_imp
                
                interp_imp_low = np.interp(rp_main, rp, impacts[model][loc][wl][1])
                interp_imp_low[rp_main>np.max(rp)] = np.nan
                med_low[:,m] = interp_imp_low
                
                interp_imp_high = np.interp(rp_main, rp, impacts[model][loc][wl][2])
                interp_imp_high[rp_main>np.max(rp)] = np.nan
                med_high[:,m] = interp_imp_high
                
            median_imp = np.nanmedian(med, axis=1)
            median_imp_low = np.nanmedian(med_low, axis=1)
            median_imp_high = np.nanmedian(med_high, axis=1)
            
            imp_w[wl] = np.vstack((median_imp, median_imp_low, median_imp_high))
        imp_loc[loc] = imp_w
    impacts['Median'] = imp_loc
    return_period['Median'] = return_period[max_rp_model]
        
    return impacts, return_period

def add_median_changes(changes):
    # data structure
    models = list(changes.keys())
    warming_levels = list(changes[models[0]].keys())
    n_cities = len(changes[models[0]][warming_levels[0]])
    imp_w = dict()
    for wl in warming_levels:
        med = np.empty([len(models), n_cities])
        for m, model in enumerate(models):
            med[m,:] = changes[model][wl]
        imp_w[wl] = np.nanmedian(med, axis=0)
    changes['Median'] = imp_w
    return changes

impacts, return_period = add_median_imp(impacts, return_period)
frac_imp, return_period = add_median_imp(frac_imp, return_period)

impacts, return_period = add_median_imp(impacts, return_period)
RP_changes = add_median_changes(RP_changes)
Pct_mort = add_median_changes(Pct_mort)
Imp_changes_rel = add_median_changes(Imp_changes_rel)
Imp_changes_pct = add_median_changes(Imp_changes_pct)


'''------------------------------------------------------------------------'''
''' save '''
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'impacts.pkl'), 'wb') as f:
    pickle.dump(impacts, f)

with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'return_period.pkl'), 'wb') as f:
    pickle.dump(return_period, f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'RP_changes.pkl'), 'wb') as f:
    pickle.dump(RP_changes, f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'Pct_mort.pkl'), 'wb') as f:
    pickle.dump(Pct_mort, f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'Imp_changes_rel.pkl'), 'wb') as f:
    pickle.dump(Imp_changes_rel, f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'Imp_changes_pct.pkl'), 'wb') as f:
    pickle.dump(Imp_changes_pct, f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'Imp_obs.pkl'), 'wb') as f:
    pickle.dump(imp_obs, f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'RP_obs.pkl'), 'wb') as f:
    pickle.dump(rp_obs, f)
    
''' impact shares '''
frac_imp = dict()
LE_models = ['CESM12-LE', 'CESM1-CAM5', 'CanESM2', 'GFDL-ESM2M', 'CSIRO-Mk3-6-0', 'Median'] # exclude GFDL-CM3 due to member size

for model in LE_models:
    print(model)
    imp_m = dict()
    
    for cc, city in enumerate(cities):
        if np.mod(cc, 25)==0: print(cc)
        mort = pd.read_excel(all_mort, sheet_name=city)
        mort_clean = _clean_data_to_365day_year(mort)
        # make annual series of mean daily mortality
        n_years = int(mort_clean.shape[0]/365)
        mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
        mort_doy_mean = np.nanmean(mort_doy,0)
        annual_mort = np.nansum(mort_doy_mean)
        
        imp_c = dict()
        
        for w, wl in enumerate(warming_levels):
            imp_c[str(wl)] = impacts[model][city][str(wl)]/annual_mort
        imp_m[city] = imp_c
    frac_imp[model] = imp_m

with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'impacts_fraction.pkl'), 'wb') as f:
    pickle.dump(frac_imp, f)


frac_imp_obs = dict()
for cc, city in enumerate(cities):
    if np.mod(cc, 25)==0: print(cc)
    RR = pd.read_excel(all_rr, sheet_name=city)
    mort = pd.read_excel(all_mort, sheet_name=city)
    frac_o, rp_o = calc_heat_mortality_fraction_obs(mort, RR,
                                all_tmm.TMM.values[all_tmm.city==city][0])
    frac_imp_obs[city] = frac_o

with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'impacts_fraction_obs.pkl'), 'wb') as f:
    pickle.dump(frac_imp_obs, f)

