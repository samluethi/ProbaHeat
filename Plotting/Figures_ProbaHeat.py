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

# all relevant plotting functions are in Plot_functions.py

'''------------------------------------------------------------------------'''
''' load all relevant data '''
DP = Path('/Users/sam/OneDrive - ETH Zurich/WCR/Projects/Heat/Data/')

with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'impacts.pkl'), 'rb') as f:
    impacts = pickle.load(f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'impacts_fraction.pkl'), 'rb') as f:
    frac_imp = pickle.load(f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'return_period.pkl'), 'rb') as f:
    return_period = pickle.load(f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'RP_changes.pkl'), 'rb') as f:
    RP_changes = pickle.load(f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'Pct_mort.pkl'), 'rb') as f:
    Pct_mort = pickle.load(f)

with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'Imp_obs.pkl'), 'rb') as f:
    imp_obs = pickle.load(f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'RP_obs.pkl'), 'rb') as f:
    rp_obs = pickle.load(f)
with open(os.path.join(DP, 'ESM/Global/run_220322/Results/', 'impacts_fraction_obs.pkl'), 'rb') as f:
    frac_imp_obs = pickle.load(f)


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
LE_models = ['CESM12-LE', 'CESM1-CAM5', 'CanESM2', 'GFDL-ESM2M', 'CSIRO-Mk3-6-0'] # exclude GFDL-CM3 due to member size

'''------------------------------------------------------------------------'''
''' clean for uncertain location '''
def set_uncertainty_exlude(impacts, cities=cities, LE_models=LE_models):

    unc = np.zeros([len(cities), len(LE_models)])
    for m, model in enumerate(LE_models):
        rep = return_period[model]
        for l, loc in enumerate(cities):
            imp = impacts[model][loc][str('0.7')]

            fit = np.interp(100., rep, imp[0,:])
            low = np.interp(100., rep, imp[1,:])
            high = np.interp(100., rep, imp[2,:])
            
            if ((high-low) > 10*fit) & (fit!=high!=low):
                unc[l,m] = 1
    unc_ex = np.sum(unc,1)==0
    
    return unc_ex

unc_ex = set_uncertainty_exlude(impacts)


'''------------------------------------------------------------------------'''
''' Draft plots'''
'''------------------------------------------------------------------------'''
''' Figure 1 - Grid plot '''
locloc = ['spal.bra9718', 'pars.fra0015', 'bngk.tha9908']

loc_names = []
for loc in locloc:
    loc_names.append(city_meta.cityname.values[city_meta.citiescode==loc][0])

obs_range_loc = []
for loc in locloc:
    mort = pd.read_excel(all_mort, sheet_name=loc)
    obs_range_loc.append([mort.date.min().year, mort.date.max().year])
    
imp_rank_loc = []
for loc in locloc:
    mort = pd.read_excel(all_mort, sheet_name=loc)
    RR = pd.read_excel(all_rr, sheet_name=loc)
    TMM = all_tmm.TMM.values[all_tmm.city==loc][0]
    rank_imp, rank_year = get_rank_heat_mort_obs(mort, RR, TMM)
    imp_rank_loc.append(rank_year[-1])

plot_grid_figure(locloc, loc_names, obs_range_loc, all_rr, all_mort,
                      frac_imp, return_period, frac_imp_obs, rp_obs, frac=True)

'''------------------------------------------------------------------------'''
''' Figure 2 - Return period map '''

# make colormap
interval = np.linspace(0.1, 1.0)
colors = plt.cm.magma(interval)
cmap = LinearSegmentedColormap.from_list('name', colors)

plot_rp_change_map_log(RP_changes, city_meta,
                   model='Median', wl=1.2, cmap=cmap, bar=False, unc_ex=unc_ex)
plot_rp_change_map_log(RP_changes, city_meta,
                   model='Median', wl=1.5, cmap=cmap, bar=False, unc_ex=unc_ex)
plot_rp_change_map_log(RP_changes, city_meta,
                   model='Median', wl=2.0, cmap=cmap, bar=False, unc_ex=unc_ex)
plot_rp_change_map_log(RP_changes, city_meta,
                   model='Median', wl=2.0, cmap=cmap, bar=True, barpos='low')
plot_rp_change_map_log(RP_changes, city_meta,
                   model='Median', wl=2.0, cmap=cmap, bar=True, barpos='right',
                   unc_ex=unc_ex)

'''------------------------------------------------------------------------'''
''' Figure 3 - Impact of extremes map '''
# make colormap
interval = np.linspace(0.2, 1.0)
colors = plt.cm.cividis(interval)
cmap = LinearSegmentedColormap.from_list('name', colors)


ticks=[1,2,3,4,5,7,10,15]
plot_mort_change_map_log(Pct_mort, city_meta,
                         title='Share of heat mortality of a 100-year season',
                         subtitle='Under 0.7$^\circ$C warming', cmax=14,
                         model='Median', wl=0.7,
                         scale_label="Share of heat mortality (%)",
                         cmap='cividis', bar=False, barpos='right', ticks=ticks)
plot_mort_change_map_log(Pct_mort, city_meta,
                         title='Share of heat mortality of a 100-year season',
                         subtitle='Under 1.2$^\circ$C warming', cmax=14,
                         model='Median', wl=1.2,
                         scale_label="Share of heat mortality (%)",
                         cmap='cividis', bar=False, barpos='right', ticks=ticks)
plot_mort_change_map_log(Pct_mort, city_meta,
                         title='Share of heat mortality of a 100-year season',
                         subtitle='Under 1.5$^\circ$C warming', cmax=14,
                         model='Median', wl=1.5,
                         scale_label="Share of heat mortality (%)",
                         cmap='cividis', bar=False, barpos='right', ticks=ticks)
plot_mort_change_map_log(Pct_mort, city_meta,
                         title='Share of heat mortality of a 100-year season',
                         subtitle='Under 2.0$^\circ$C warming', cmax=14,
                         model='Median', wl=2.0,
                         scale_label="Share of heat mortality (%)",
                         cmap='cividis', bar=False, barpos='right', ticks=ticks)
plot_mort_change_map_log(Pct_mort, city_meta,
                         title='Share of heat mortality of a 100-year season',
                         subtitle='Under 2.0$^\circ$C warming', cmax=14,
                         model='Median', wl=2.0,
                         scale_label="Fraction of heat mortality (%)",
                         cmap='cividis', bar=True, barpos='low', ticks=ticks)

'''------------------------------------------------------------------------'''
''' Figure 4 - Petri dish '''
plot_petri_dish(impacts, return_period, location='pars.fra0015',
                model='Median', title='', sort_points=False)

'''------------------------------------------------------------------------'''
''' Appendix plot '''

''' Uncertainty of discrete return periods per model '''
plot_grid_figure_WL_RP(locloc, loc_names, impacts, return_period, ylim=[5000, 15000, 2500])

''' Stochastic uncertainty of ExFreqCurve '''
dat_path= DP.joinpath('ESM/Global/run_220322/')
plot_grid_line_ref_bootstrap(dat_path, LE_models, locloc, imp_obs, rp_obs,
                             all_rr, all_mort, all_tmm)

''' Impact exceedance frequency curve per model '''
plot_all_model_grid(locloc, impacts, return_period,
                    loc_title=loc_names, ylim=[5000, 15000, 2500])

''' Changes in RP & extreme impacts per model '''
plot_grid_figure_model_change(locloc, locloc, RP_changes, Pct_mort,
                              ylim_RP=[60,50,30], ylim_pct=[6,24,10])

'''------------------------------------------------------------------------'''
''' Table for Supplementary '''

# get locs IQR
RP_changes_IQR = dict()
pct_mort_IQR = dict()
for w, wl in enumerate(warming_levels):
    rp_wl = np.zeros([len(cities),len(LE_models)])
    pm_wl = np.zeros([len(cities),len(LE_models)])
    rpiqr = dict()
    pmiqr = dict()
    for m, model in enumerate(LE_models):
        rp_wl[:,m] = RP_changes[model][wl]
        pm_wl[:,m] = Pct_mort[model][wl]
        
    rpiqr['lower'] = np.sort(rp_wl)[:,-4:-3]
    rpiqr['upper'] = np.sort(rp_wl)[:,-2:-1]
    pmiqr['lower'] = np.sort(pm_wl)[:,-4:-3]
    pmiqr['upper'] = np.sort(pm_wl)[:,-2:-1]

    RP_changes_IQR[wl] = rpiqr
    pct_mort_IQR[wl] = pmiqr
    
''' make table '''

results_all_loc = pd.DataFrame()

results_all_loc['Location'] = city_meta.cityname
results_all_loc['Country'] = city_meta.countryname
results_all_loc['RP_2020'] = RP_changes['Median'][1.2]
results_all_loc['RP_2020_low'] = RP_changes_IQR[1.2]['lower']
results_all_loc['RP_2020_high'] = RP_changes_IQR[1.2]['upper']
results_all_loc['RP_15'] = RP_changes['Median'][1.5]
results_all_loc['RP_15_low'] = RP_changes_IQR[1.5]['lower']
results_all_loc['RP_15_high'] = RP_changes_IQR[1.5]['upper']
results_all_loc['RP_20'] = RP_changes['Median'][2.0]
results_all_loc['RP_20_low'] = RP_changes_IQR[2.0]['lower']
results_all_loc['RP_20_high'] = RP_changes_IQR[2.0]['upper']

results_all_loc['PM_2000'] = Pct_mort['Median'][0.7]
results_all_loc['PM_2000_low'] = pct_mort_IQR[0.7]['lower']
results_all_loc['PM_2000_high'] = pct_mort_IQR[0.7]['upper']
results_all_loc['PM_2020'] = Pct_mort['Median'][1.2]
results_all_loc['PM_2020_low'] = pct_mort_IQR[1.2]['lower']
results_all_loc['PM_2020_high'] = pct_mort_IQR[1.2]['upper']
results_all_loc['PM_15'] = Pct_mort['Median'][1.5]
results_all_loc['PM_15_low'] = pct_mort_IQR[1.5]['lower']
results_all_loc['PM_15_high'] = pct_mort_IQR[1.5]['upper']
results_all_loc['PM_20'] = Pct_mort['Median'][2.0]
results_all_loc['PM_20_low'] = pct_mort_IQR[2.0]['lower']
results_all_loc['PM_20_high'] = pct_mort_IQR[2.0]['upper']

results_all_loc.to_excel(os.path.join(DP, 'Results_All_Locations.xlsx'))
