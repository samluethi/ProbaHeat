#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 17 2022
@author: Samuel Luethi - samuel.luethi@usys.ethz.ch

"""
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.image as mpimg
from matplotlib.patches import Circle
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LogNorm


import cartopy.crs as crs
import cartopy.feature as cfeature
import seaborn as sns
import numpy as np
import pandas as pd

# This script contains all plotting functions.
# See script Figures_ProbaHeat.py for running the functions

'''------------------------------------------------------------------------'''
''' Grid plot Figure 1 '''
def plot_grid_figure(locations, loc_names, obs_range_loc, all_rr, all_mort,
                     impacts, return_period, imp_obs, rp_obs,
                     wl=[0.7, 1.2, 1.5, 2.0], frac=False):

    plt.figure(figsize=(20, 20))
    grid = plt.GridSpec(3, len(locations), wspace=0.3, hspace=0.6) # create a grid for the subplots
    
    for l, loc in enumerate(locations):
        
        if l==0: ylab=True
        else: ylab=False
        
        ax = plt.subplot(grid[0,l])
        RR = pd.read_excel(all_rr, sheet_name=loc)
        mort = pd.read_excel(all_mort, sheet_name=loc)
        plot_RR(RR, mort.temp, loc_names[l], col='slategray', ylab=ylab, labelfont=24)
        # legends
        if l==1:
            plt.plot(np.nan, np.nan)
            p1 = plt.plot(np.NaN, np.NaN, color="slategray", linewidth=4)
            p2 = plt.plot(np.NaN, np.NaN, '--', color="slategray", linewidth=4)
            p3 = plt.bar(np.NaN, np.NaN, color="slategray", alpha=0.8)
            p4 = plt.plot(np.NaN, np.NaN, linestyle='dotted',color="black")
        
            plt.legend([p1[0], p2[0],p3[0],p4[0]],
                       ['Relative risk',
                        'Extrapolated range',
                        '95% uncertainty interval',
                        'p99 of temperature'],
                       bbox_to_anchor=(0.5, -0.5), loc='lower center', ncol=4,
                       fontsize=24, frameon=False)
                
        ax = plt.subplot(grid[1,l])
        plot_line_ex_freq(impacts, return_period, imp_obs, rp_obs, loc,
                                   title='', fill=False, cmap='cividis', model='Median',
                                   wl=[0.7, 1.2, 1.5, 2.0], xmax=500, titlefont=32,
                                   tickfont=24, labelfont=24, ylab=ylab,
                                   obs_leg=True, obs_range=obs_range_loc[l],
                                   err_bar=True, frac=frac)
        # legends
        if l==1:
            cm = plt.get_cmap('cividis')
            cols = [cm(i / len(wl)) for i in range(len(wl))]
            p1 = plt.plot(np.NaN, np.NaN, marker='d', color="black", markersize=12, linewidth=0)
            p2 = plt.plot(np.NaN, np.NaN, color=cols[0], linewidth=4)
            p3 = plt.plot(np.NaN, np.NaN, color=cols[1], linewidth=4)
            p4 = plt.plot(np.NaN, np.NaN, color=cols[2], linewidth=4)
            p5 = plt.plot(np.NaN, np.NaN, color=cols[3], linewidth=4)
    
            plt.legend([p1[0], p2[0],p3[0], p4[0], p5[0]],
                       ['Obs. estimate',
                        'Risk 2000',
                        'Risk 2020',
                        'Risk at 1.5$^\circ$C',
                        'Risk at 2.0$^\circ$C'],
                       bbox_to_anchor=(0.5, -0.5), loc='lower center', ncol=5,
                       fontsize=24, frameon=False)
        
        ax = plt.subplot(grid[2,l])
        plot_dat_bar(loc, impacts, return_period, cmap='twilight', ylab=ylab,
                     labelfont=24, frac=frac)
        # legends
        if l==1:
            cm = plt.get_cmap('twilight')
            cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
            p1 = plt.bar(np.NaN, np.NaN, color=cols[1])
            p2 = plt.bar(np.NaN, np.NaN, color=cols[2])
            p3 = plt.bar(np.NaN, np.NaN, color=cols[3])
            p4 = plt.bar(np.NaN, np.NaN, color=cols[4])
            p5 = plt.bar(np.NaN, np.NaN, color=cols[5])
    
            plt.legend([p1[0], p2[0],p3[0], p4[0], p5[0]],
                       LE_models,
                       bbox_to_anchor=(0.5, -0.5), loc='lower center', ncol=5,
                       fontsize=24, frameon=False)
        


'''------------------------------------------------------------------------'''
''' Figure 2 - Return period map '''

def plot_rp_change_map(change, city_info, title='', subtitle='', cmax=50,
                       model='Median', wl=1.5, cmap='magma_r', bar=True,
                       barpos='right', unc_ex=0):
    rp_change = change[model][wl]
    lat = city_info.lat
    lon = city_info.long
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(rp_change):
        n_lat = lat[~unc_ex]
        n_lon = lon[~unc_ex]
        
        rp_change = rp_change[unc_ex]
        lat = lat[unc_ex]
        lon = lon[unc_ex]
        
    fig = plt.figure(figsize=(32,24))
    plt.rcParams.update({'font.family':'arial'})
    
    ax = fig.add_subplot(1,1,1, projection=crs.Robinson())
    
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND, color="lightgrey", alpha=0.5)
    ax.add_feature(cfeature.OCEAN, color="skyblue", alpha=0.4)
    ax.gridlines()
    
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(rp_change):
        plt.scatter(x=n_lon, y=n_lat, c='grey', cmap=cmap,
                    marker=".",
                    s=600,
                    alpha=0.5,
                    transform=crs.PlateCarree())
        
    sc = plt.scatter(x=lon, y=lat, c=rp_change, cmap=cmap,
                marker=".",
                s=600,
                alpha=1.0,
                transform=crs.PlateCarree())
    
    ax.set_extent([180, -180, -90, 90], crs.PlateCarree())
    
    if bar:
        if barpos=='right':
            cbax = make_axes_locatable(ax).append_axes(
                'right', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='vertical', extend='neither')
        elif barpos=='low':
            cbax = make_axes_locatable(ax).append_axes(
                'bottom', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='horizontal', extend='neither')
        cbar.set_label("New return period (Year)", fontsize=48)
        cbar.ax.tick_params(labelsize=42)
    plt.clim(0, cmax)

    ax.set_title(f"{title}\n{subtitle}",fontsize=36)
    plt.show()
    
def plot_rp_change_map_log(change, city_info, title='', subtitle='', cmax=50,
                       model='Median', wl=1.5, cmap='magma_r', bar=True,
                       barpos='right', unc_ex=0, ticks=[1,2,5,10,20,50]):
    rp_change = change[model][wl]
    lat = city_info.lat
    lon = city_info.long
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(rp_change):
        n_lat = lat[~unc_ex]
        n_lon = lon[~unc_ex]
        
        rp_change = rp_change[unc_ex]
        lat = lat[unc_ex]
        lon = lon[unc_ex]
        
    fig = plt.figure(figsize=(32,24))
    plt.rcParams.update({'font.family':'arial'})
    
    ax = fig.add_subplot(1,1,1, projection=crs.Robinson())
    
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND, color="lightgrey", alpha=0.5)
    ax.add_feature(cfeature.OCEAN, color="skyblue", alpha=0.4)
    ax.gridlines()
    
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(rp_change):
        plt.scatter(x=n_lon, y=n_lat, c='grey',
                    marker=".",
                    s=600,
                    alpha=0.5,
                    transform=crs.PlateCarree())
        
    sc = plt.scatter(x=lon, y=lat, c=rp_change, cmap=cmap,
                norm=LogNorm(vmin=ticks[0], vmax=ticks[-1]),
                marker=".",
                s=600,
                alpha=1.0,
                transform=crs.PlateCarree())
    
    ax.set_extent([180, -180, -90, 90], crs.PlateCarree())
    
    if bar:
        if barpos=='right':
            cbax = make_axes_locatable(ax).append_axes(
                'right', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='vertical', extend='neither')
            cbar.set_ticks(ticks)
            cbar.ax.set_yticklabels(ticks)
        elif barpos=='low':
            cbax = make_axes_locatable(ax).append_axes(
                'bottom', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='horizontal', extend='neither')
            cbar.set_ticks(ticks)
            cbar.ax.set_xticklabels(ticks)
        cbar.set_label("New return period (Year)", fontsize=48)
        cbar.ax.tick_params(labelsize=42)
    plt.clim(ticks[0], ticks[-1])

    ax.set_title(f"{title}\n{subtitle}",fontsize=36)
    plt.show()

def plot_mort_change_map(change, city_info, title='', subtitle='', cmax=250,
                       model='Median', wl=1.5, scale_label="Increase in mortality [%]",
                       cmap='viridis',bar=True, barpos='right', unc_ex=0):
    
    mort_change = change[model][wl]
    lat = city_info.lat
    lon = city_info.long
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(change):
        n_lat = lat[~unc_ex]
        n_lon = lon[~unc_ex]
        
        mort_change = mort_change[unc_ex]
        lat = lat[unc_ex]
        lon = lon[unc_ex]
        
    fig = plt.figure(figsize=(32,24))
    plt.rcParams.update({'font.family':'arial'})
    
    ax = fig.add_subplot(1,1,1, projection=crs.Robinson())
    
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND, color="lightgrey", alpha=0.5)
    ax.add_feature(cfeature.OCEAN, color="skyblue", alpha=0.4)
    ax.gridlines()
    
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(change):
        plt.scatter(x=n_lon, y=n_lat, c='grey', cmap=cmap,
                    marker=".",
                    s=600,
                    alpha=0.25,
                    transform=crs.PlateCarree())
    sc = plt.scatter(x=lon, y=lat, c=mort_change*100, cmap=cmap,
                marker=".",
                s=600,
                alpha=1.0,
                transform=crs.PlateCarree())
    
    ax.set_extent([180, -180, -90, 90], crs.PlateCarree())
    
    if bar:
        if barpos=='right':
            cbax = make_axes_locatable(ax).append_axes(
                'right', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='vertical', extend='neither')
        elif barpos=='low':
            cbax = make_axes_locatable(ax).append_axes(
                'bottom', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='horizontal', extend='neither')
        cbar.set_label(scale_label, fontsize=48)
        cbar.ax.tick_params(labelsize=42)
    plt.clim(0, cmax)
    
    plt.show()

def plot_mort_change_map_log(change, city_info, title='', subtitle='', cmax=250,
                       model='Median', wl=1.5, scale_label="Increase in mortality [%]",
                       cmap='viridis',bar=True, barpos='right', unc_ex=0, ticks=[1,2,5,10,20,50]):
    
    mort_change = change[model][wl]
    lat = city_info.lat
    lon = city_info.long
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(change):
        n_lat = lat[~unc_ex]
        n_lon = lon[~unc_ex]
        
        mort_change = mort_change[unc_ex]
        lat = lat[unc_ex]
        lon = lon[unc_ex]
        
    fig = plt.figure(figsize=(32,24))
    plt.rcParams.update({'font.family':'arial'})
    
    ax = fig.add_subplot(1,1,1, projection=crs.Robinson())
    
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAND, color="lightgrey", alpha=0.5)
    ax.add_feature(cfeature.OCEAN, color="skyblue", alpha=0.4)
    ax.gridlines()
    
    if np.sum(unc_ex) > 0 & np.sum(unc_ex)!=len(change):
        plt.scatter(x=n_lon, y=n_lat, c='grey', cmap=cmap,
                    marker=".",
                    s=600,
                    alpha=0.25,
                    transform=crs.PlateCarree())
    sc = plt.scatter(x=lon, y=lat, c=mort_change*100, cmap=cmap,
                norm=LogNorm(vmin=ticks[0], vmax=ticks[-1]),
                marker=".",
                s=600,
                alpha=1.0,
                transform=crs.PlateCarree())
    
    ax.set_extent([180, -180, -90, 90], crs.PlateCarree())
    
    if bar:
        if barpos=='right':
            cbax = make_axes_locatable(ax).append_axes(
                'right', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='vertical', extend='neither')
            cbar.set_ticks(ticks)
            cbar.ax.set_yticklabels(ticks)
        elif barpos=='low':
            cbax = make_axes_locatable(ax).append_axes(
                'bottom', size="6.5%", pad=0.1, axes_class=plt.Axes)
            cbar = plt.colorbar(sc, cax=cbax, orientation='horizontal', extend='neither')
            cbar.set_ticks(ticks)
            cbar.ax.set_xticklabels(ticks)
        cbar.set_label(scale_label, fontsize=48)
        cbar.ax.tick_params(labelsize=42)
    plt.clim(ticks[0], ticks[-1])
    
    plt.show()

'''------------------------------------------------------------------------'''
''' line plot '''
def plot_line_ex_freq(imp_mod, rp_mod, imp_obs, rp_obs, city,
                           title='', fill=False, cmap='cividis', model='Median',
                           wl=[0.7, 1.2, 1.5, 2.0], xmax=500, titlefont=32,
                           tickfont=20, labelfont=22, ylab=True,
                           obs_leg=False, obs_range=[np.nan, np.nan],
                           err_bar=False, RP=np.nan, frac=False,
                           obs_year=False, obs_year_max=np.nan):
    cm = plt.get_cmap(cmap)
    cols = [cm(i / len(wl)) for i in range(len(wl))]
    
    if frac:
        for w, wlevel in enumerate(wl):
            plt.plot(rp_mod[model], imp_mod[model][city][str(wlevel)][0]*100, color=cols[w], linewidth=3)
        plt.plot(rp_obs[city], imp_obs[city][0]*100, color='black', linewidth=2.0)
        plt.scatter(rp_obs[city][rp_obs[city]>2], imp_obs[city][0][rp_obs[city]>2]*100, marker='d', color='black', s=80)    
        ymax = np.interp(xmax, rp_mod[model],imp_mod[model][city][str(wlevel)][0])*1.1*100
    else:
        for w, wlevel in enumerate(wl):
            plt.plot(rp_mod[model], imp_mod[model][city][str(wlevel)][0], color=cols[w], linewidth=3)
        plt.plot(rp_obs[city], imp_obs[city][0], color='black', linewidth=2.0)
        plt.scatter(rp_obs[city][rp_obs[city]>2], imp_obs[city][0][rp_obs[city]>2], marker='d', color='black', s=80)
        ymax = np.interp(xmax, rp_mod[model],imp_mod[model][city][str(wlevel)][0])*1.1
    
    plt.title(title, fontsize=titlefont)
    plt.xscale('log')
    plt.xticks(fontsize=tickfont)
    plt.yticks(fontsize=tickfont)
    plt.xlim(1,xmax)
    plt.ylim(0,ymax)
    sns.despine()
    plt.xlabel("Return period (Year)", fontsize=labelfont)
    if ylab:
        if frac:
            plt.ylabel("Heat mortality (%)", fontsize=labelfont)
        else:
            plt.ylabel("Heat mortality (#)", fontsize=labelfont)
    if err_bar:
        for w, wlevel in enumerate(wl):
            imp_fit = np.interp(RP, rp_mod[model], imp_mod[model][city][str(wlevel)][0])
            imp_low = np.interp(RP, rp_mod[model], imp_mod[model][city][str(wlevel)][1])
            imp_high = np.interp(RP, rp_mod[model], imp_mod[model][city][str(wlevel)][2])
            err=np.array(imp_fit-imp_low, imp_high-imp_fit)
            plt.errorbar(xmax+50+w*100, imp_fit, yerr=err, fmt='+', color=cols[w],
                         clip_on=False,  markersize='10', capsize=1, elinewidth=3)
        
    if obs_leg:
        p1 = plt.plot(np.NaN, np.NaN, color="black", linewidth=4)
        lege = plt.legend([p1[0]],
               [str(obs_range[0])+'-'+str(obs_range[1])], loc='upper left',
                fontsize=20, frameon=False)
        plt.gca().add_artist(lege)
    
    if obs_year:
        plt.text(rp_obs[city][-1]*0.9, imp_obs[city][0][-1]*0.8, str(obs_year_max),
                 fontsize=16)

#plot_line_ex_freq(impacts, return_period, imp_obs, rp_obs, locations[10])    
  

# plot_line_all_to_grid(plot_line_all_models, impacts_w12, return_period, imp_obs, rp_obs,
#                       locations, loc_names, y_max_loc, title='', obs_leg=True, obs_range=obs_range_loc)

def plot_line_all_to_grid(plot_fun, imp_mod, rp_mod, imp_obs, rp_obs, locations,
                          loc_title, ymaxloc, title='', wl=1.2, obs_leg=False, obs_range=[np.nan, np.nan],
                          cmap='viridis', **kwargs):

    
    plt.figure(figsize=(24, 20))
    plt.suptitle(title, fontsize=44)

    
    for i in range(12):
        plt.subplot(5,3,i+1)
        plt.subplots_adjust(hspace=.8)
        plt.subplots_adjust(wspace=.3)
        plot_fun(imp_mod, rp_mod, imp_obs, rp_obs, locations[i],
                 title=loc_title[i], wl=wl, obs_leg=obs_leg, obs_range=obs_range[i],
                 cmap=cmap, ymax=ymaxloc[i])
    
    cm = plt.get_cmap(cmap)
    cols = [cm(i / len(imp_mod)) for i in range(len(imp_mod))]
    models = imp_mod.keys()
    pp = []
    pp.append(plt.plot(np.NaN, np.NaN, color="black", linewidth=4)[0])
    pp.append(plt.plot(np.NaN, np.NaN, color="black", linestyle='dashed', linewidth=4)[0])
    for m in range(len(models)-1):
        pp.append(plt.plot(np.NaN, np.NaN, color=cols[m], linewidth=4)[0])
    leg = []
    leg.append('Observational estimate')
    leg.append('Median')
    for model in models:
        if not model=='Median':
            leg.append(model)
    plt.legend(pp, leg,
               bbox_to_anchor=(-1.0, -1.5), loc='lower center', ncol=3,
               fontsize=24, frameon=False)
    
# plot_line_all_to_grid(plot_line_all_models, impacts, return_period, imp_obs, rp_obs,
#                       locations, loc_names, y_max_loc, title='', wl=1.5, obs_leg=True, obs_range=obs_range_loc)
'''------------------------------------------------------------------------'''
''' error bar plots '''
def plot_dat_bar(loc, impacts, return_period, title='', cmap='viridis',
                      LE_models=LE_models, labelfont=22, titlefont=22,
                      ylab=True, x_label='', frac=False, short_name=False):
    
    cm = plt.get_cmap(cmap)
    cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
    if short_name:
        warming=np.array(['00', '20', '1.5', '2.0'])
    else:
        warming=np.array([str(2000), '2020', '1.5$^\circ$C', '2.0$^\circ$C'])
    
    RP_imp = prep_dat_bar(loc, impacts, return_period)
    x = np.arange(0, len(warming))
    if frac:
        for m, model in enumerate(LE_models):
            width = 0.1
            x_fix = x - 2*width + m*width
            err = [np.array(RP_imp['fit'][model])*100-np.array(RP_imp['low'][model])*100,
                         np.array(RP_imp['high'][model])*100-np.array(RP_imp['fit'][model])*100]
            plt.bar(x_fix, np.array(RP_imp['fit'][model])*100, width, color=cols[m+1], yerr=err, ecolor='gray')
        plt.errorbar(x, np.array(RP_imp['fit']['Median'])*100, fmt='_', color='black',
                     markersize='60', markeredgewidth=2, elinewidth=3)
    else:
        for m, model in enumerate(LE_models):
            width = 0.1
            x_fix = x - 2*width + m*width
            err = [np.array(RP_imp['fit'][model])-np.array(RP_imp['low'][model]),
                         np.array(RP_imp['high'][model])-np.array(RP_imp['fit'][model])]
            plt.bar(x_fix, RP_imp['fit'][model], width, color=cols[m+1], yerr=err, ecolor='gray')
        plt.errorbar(x, RP_imp['fit']['Median'], fmt='_', color='black',
                     markersize='60', markeredgewidth=2, elinewidth=3)
        
    if ylab:
        if frac:
            plt.ylabel('Heat Mortality (%)', fontsize=labelfont)
        else:
            plt.ylabel('Heat Mortality (#)', fontsize=labelfont)
    plt.xticks(x, warming, fontsize=labelfont)
    plt.xlabel(x_label, fontsize=labelfont)
    plt.yticks(fontsize=labelfont)
    sns.despine()
    plt.title(title, fontsize=titlefont)



def prep_dat_bar(loc, impacts, return_period,
                 wl=np.array([0.7, 1.2, 1.5, 2.0]), RP=100.):
    RP_imp = dict()
    RP_fit = dict()
    RP_low = dict()
    RP_high = dict()
    
    for model, m_imp in impacts.items():
        fit = list()
        low = list()
        high = list()
        for w, wlevel in enumerate(wl):
            rp_mod = return_period[model]
            imp = m_imp[loc][str(wlevel)]
            
            fit.append(np.interp(RP, rp_mod, imp[0,:]))
            low.append(np.interp(RP, rp_mod, imp[1,:]))
            high.append(np.interp(RP, rp_mod, imp[2,:]))
        RP_fit[model] = fit
        RP_low[model] = low
        RP_high[model] = high
    RP_imp = dict()
    RP_imp['fit'] = RP_fit
    RP_imp['low'] = RP_low
    RP_imp['high'] = RP_high
    
    return RP_imp

# plot_dat_bar(loc, impacts, return_period, cmap='jet')
# plot_dat_errorbar(loc, impacts, return_period, cmap='jet')

'''------------------------------------------------------------------------'''
''' RR plots '''
def plot_RR(RR, T, city_name='', col='indianred', col2='black',
            include_HM=False, ylab=True, labelfont=22):
    
    perc_1 = T.quantile(0.01)
    perc_99 = T.quantile(0.99)
    xmin = perc_1
    xmax = T.max()+2.
    ymin = 0.9
    ymax = RR.RRfit[np.argmin(abs(RR.temp-T.max()))]*1.5
    ymax = 2.0
    
    RR_99 = RR[RR.temp.between(perc_1, T.max())]
    
    # including extrapolation
    plt.plot(RR_99.temp.values, RR_99.RRfit.values, '-', color=col)
    plt.fill_between(RR_99.temp.values, RR_99.RRlow.values, RR_99.RRhigh.values,
                    color=col, alpha=0.2)
    
                    
    plt.plot(RR.temp.values, RR.RRfit.values, '--', color=col)
    plt.fill_between(RR.temp.values, RR.RRlow.values, RR.RRhigh.values,
                    color=col, alpha=0.1)
    
    plt.fill_between(RR.temp.values, RR.RRlow.values, RR.RRhigh.values,
                    color=col, alpha=0.1)
    if include_HM:
        plt.fill_between(RR.temp.values[RR.temp.values>=TMM],
                        np.repeat(RR.RRfit.values[TMM==RR.temp.values],
                                  len(RR.temp.values[RR.temp.values>=TMM])),
                        RR.RRfit.values[RR.temp.values>=TMM],
                        color=col2, alpha=0.1)
    
    # p1 = plt.plot(np.NaN, np.NaN, color=col)
    # p2 = plt.fill(np.NaN, np.NaN, color=col, alpha=0.2)
    plt.plot([xmin, xmax], [1.0, 1.0], 'k-',linewidth=1.)
    plt.plot([perc_1, perc_1], [1.0, 10.0], 'k--',linewidth=1.)
    plt.plot([perc_99, perc_99], [1.0, 10.0], 'k--',linewidth=1.)
    plt.ylim([ymin, ymax])
    plt.xlim([xmin, xmax])   
    plt.xticks(fontsize=labelfont)
    plt.yticks(fontsize=labelfont)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    plt.xlabel("Temperature (C)", fontsize=labelfont)
    if ylab:
        plt.ylabel("Relative risk", fontsize=labelfont)
    plt.title(city_name)

# plot_RR(RR, mort.temp, 'abc', col='indianred')


'''------------------------------------------------------------------------'''
''' Appendix Plot all models '''


def plot_line_all_models(imp_mod, rp_mod, imp_obs, rp_obs, city, title='',
                         wl=1.2, cmap='viridis', ymax=0,
                         xmax=1500, titlefont=32, tickfont=20, labelfont=22,
                         obs_leg=False, obs_range=[np.nan, np.nan]):
    cm = plt.get_cmap(cmap)
    cols = [cm(i / len(imp_mod)) for i in range(len(imp_mod))]
    models = imp_mod.keys()

    for m, model in enumerate(models):
        if model=='Median':
            plt.plot(rp_mod[model], imp_mod[model][city][str(wl)][0], color='black', linestyle='dashed', linewidth=3)
        else:
            plt.plot(rp_mod[model], imp_mod[model][city][str(wl)][0], color=cols[m], linewidth=2)
    plt.plot(rp_obs[city], imp_obs[city][0], color='black', linewidth=3)

    plt.title(title, fontsize=titlefont)
    plt.xscale('log')
    plt.xticks(fontsize=tickfont)
    plt.yticks(fontsize=tickfont)
    plt.xlim(1,xmax)
    if ymax>0:
        plt.ylim(0,ymax)
    plt.xlabel("Exceedance frequency [Year]", fontsize=labelfont)
    plt.ylabel("Excess mortality [#]", fontsize=labelfont)
    if obs_leg:
        p1 = plt.plot(np.NaN, np.NaN, color="black", linewidth=4)
        plt.legend([p1[0]],
               [str(obs_range[0])+'-'+str(obs_range[1])], loc='upper left',
                fontsize=20, frameon=False)

# plot_line_all_models(impacts, return_period, imp_obs, rp_obs, locations[10])  
def plot_line_all_models_no_obs(imp_mod, rp_mod, city, title='',
                         wl=1.2, cmap='viridis', models='', ymax=0,
                         xmax=1500, titlefont=32, tickfont=20, labelfont=22,
                         ylab=True, obs_leg=False, obs_range=[np.nan, np.nan]):
    cm = plt.get_cmap(cmap)
    cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
    if models=='':
        models = imp_mod.keys()

    for m, model in enumerate(models):
        if model=='Median':
            plt.plot(rp_mod[model], imp_mod[model][city][str(wl)][0], color='black', linestyle='dashed', linewidth=3)
        else:
            plt.plot(rp_mod[model], imp_mod[model][city][str(wl)][0], color=cols[m+1], linewidth=2)

    plt.title(title, fontsize=titlefont)
    plt.xscale('log')
    plt.xticks(fontsize=tickfont)
    plt.yticks(fontsize=tickfont)
    plt.xlim(1,xmax)
    sns.despine()
    if ymax>0:
        plt.ylim(0,ymax)
    plt.xlabel("Return period (Year)", fontsize=labelfont)
    if ylab:
        plt.ylabel("Heat mortality (#)", fontsize=labelfont)
    if obs_leg:
        p1 = plt.plot(np.NaN, np.NaN, color="black", linewidth=4)
        plt.legend([p1[0]],
               [str(obs_range[0])+'-'+str(obs_range[1])], loc='upper left',
                fontsize=20, frameon=False)
        
def plot_all_model_grid(locations, imp_mod, rp_mod,
                          loc_title, title='', LE_models=LE_models, ylim=np.nan,
                          obs_leg=False, obs_range=[np.nan, np.nan],
                          cmap='twilight', model='Median', wl=[0.7, 1.2, 1.5, 2.0],
                          **kwargs):

    plt.figure(figsize=(len(locations)*6, 20))
    grid = plt.GridSpec(4, len(locations), wspace=0.5, hspace=0.4) # create a grid for the subplots
    
    for l, loc in enumerate(locations):
        if l==0: ylab=True
        else: ylab=False
        for w, warming in enumerate(wl):
            if w==0: tit=loc_title[l]
            else: tit=''

            ax = plt.subplot(grid[w,l])
            plot_line_all_models_no_obs(imp_mod, rp_mod, loc, title=tit,
                                         wl=warming, cmap=cmap, models=LE_models,
                                         ymax=ylim[l], xmax=1680, #1500
                                         titlefont=32, tickfont=20,
                                         labelfont=22, ylab=ylab)
           
        
    cm = plt.get_cmap(cmap)
    cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
    p1 = plt.plot(np.NaN, np.NaN, color=cols[1], linewidth=4)
    p2 = plt.plot(np.NaN, np.NaN, color=cols[2], linewidth=4)
    p3 = plt.plot(np.NaN, np.NaN, color=cols[3], linewidth=4)
    p4 = plt.plot(np.NaN, np.NaN, color=cols[4], linewidth=4)
    p5 = plt.plot(np.NaN, np.NaN, color=cols[5], linewidth=4)

    plt.legend([p1[0], p2[0],p3[0], p4[0], p5[0]],
               LE_models,
               bbox_to_anchor=(-1.1, -0.7), loc='lower center', ncol=5,
               fontsize=22, frameon=False)
    
# plot_all_model_grid(locloc, impacts, return_period,
#                     loc_title=loc_names, ylim=[5000, 15000, 2500])

    
'''------------------------------------------------------------------------'''
''' Appendix Plot uncertainty different warming levels & RPs'''
def plot_grid_figure_WL_RP(locations, loc_names, impacts, return_period,
                           wl=[0.7, 1.2, 1.5, 2.0], ylim=np.nan):

    plt.figure(figsize=(20, 20))
    grid = plt.GridSpec(len(wl), len(locations), wspace=0.3, hspace=0.6) # create a grid for the subplots
    
    for l, loc in enumerate(locations):
        for w, warm in enumerate(wl):
            if l==0: ylab=True
            else: ylab=False
            
            if not any(np.isnan(ylim)): ymax=ylim[l]
            else: ymax = 0
                
            
            ax = plt.subplot(grid[w,l])
            if w==0:
                plot_dat_bar_RP(loc, impacts, return_period, wl=warm, title=loc_names[l],
                                cmap='twilight', ylab=ylab, ymax=ymax, labelfont=24)
            else:
                plot_dat_bar_RP(loc, impacts, return_period, wl=warm,
                                cmap='twilight', ylab=ylab, ymax=ymax, labelfont=24)
            # legends
        if l==1:
            cm = plt.get_cmap('twilight')
            cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
            p1 = plt.bar(np.NaN, np.NaN, color=cols[1])
            p2 = plt.bar(np.NaN, np.NaN, color=cols[2])
            p3 = plt.bar(np.NaN, np.NaN, color=cols[3])
            p4 = plt.bar(np.NaN, np.NaN, color=cols[4])
            p5 = plt.bar(np.NaN, np.NaN, color=cols[5])
    
            plt.legend([p1[0], p2[0],p3[0], p4[0], p5[0]],
                       LE_models,
                       bbox_to_anchor=(0.5, -0.7), loc='lower center', ncol=5,
                       fontsize=24, frameon=False)
            
def plot_dat_bar_RP(loc, impacts, return_period, wl=1.2,
                 title='', cmap='viridis', LE_models=LE_models,
                 labelfont=22, titlefont=22, ylab=True, ymax=0):
    
    cm = plt.get_cmap(cmap)
    cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
    ret_per = np.array(['10', '50', '100', '500'])
    
    RP_imp = prep_dat_bar_RP(loc, impacts, return_period, wl=wl)
    x = np.arange(0, len(ret_per))
    for m, model in enumerate(LE_models):
        width = 0.1
        x_fix = x - 2*width + m*width
        err = [np.array(RP_imp['fit'][model])-np.array(RP_imp['low'][model]),
                     np.array(RP_imp['high'][model])-np.array(RP_imp['fit'][model])]
        plt.bar(x_fix, RP_imp['fit'][model], width, color=cols[m+1], yerr=err, ecolor='gray')
    plt.errorbar(x, RP_imp['fit']['Median'], fmt='_', color='black',
                 markersize='60', markeredgewidth=2, elinewidth=3)
    if ylab:
        plt.ylabel('Heat Mortality (#)', fontsize=labelfont)
    if ymax>0:
        plt.ylim(0,ymax)
    plt.xticks(x, ret_per, fontsize=labelfont)
    plt.xlabel('Return period (year)', fontsize=labelfont)
    plt.yticks(fontsize=labelfont)
    sns.despine()
    plt.title(title, fontsize=titlefont)

def prep_dat_bar_RP(loc, impacts, return_period,
                 RP=np.array([10, 50, 100, 500]), wl=1.2):
    RP_imp = dict()
    RP_fit = dict()
    RP_low = dict()
    RP_high = dict()
    
    for model, m_imp in impacts.items():
        fit = list()
        low = list()
        high = list()
        for r, rp in enumerate(RP):
            rp_mod = return_period[model]
            imp = m_imp[loc][str(wl)]
            
            fit.append(np.interp(rp, rp_mod, imp[0,:]))
            low.append(np.interp(rp, rp_mod, imp[1,:]))
            high.append(np.interp(rp, rp_mod, imp[2,:]))
        RP_fit[model] = fit
        RP_low[model] = low
        RP_high[model] = high
    RP_imp = dict()
    RP_imp['fit'] = RP_fit
    RP_imp['low'] = RP_low
    RP_imp['high'] = RP_high
    
    return RP_imp
'''------------------------------------------------------------------------'''
''' cat plot uncertainty '''
def prep_dat_cat_RP(impacts, loc, return_period, RP = np.array([10, 50, 100, 500]), wl=1.2):
    plt_dat = pd.DataFrame(columns=['Model', 'Return_period', 'RR', 'Value'])
    for model, m_imp in impacts.items():
        imp = m_imp[loc][str(wl)]
        rp_mod = return_period[model]
        imp_fit = np.interp(RP, rp_mod, imp[0,:])
        imp_low = np.interp(RP, rp_mod, imp[1,:])
        imp_high = np.interp(RP, rp_mod, imp[2,:])
        
        mod_dat = pd.DataFrame(columns=['Model', 'Return_period', 'RR', 'Value'])
        mod_dat['Model'] = np.repeat(model, len(RP)*3)
        mod_dat['Return_period'] = np.tile(RP, 3)
        mod_dat['RR'] = np.repeat(np.array(['fit', 'low', 'high']), len(RP))
        mod_dat['Value'] = np.ravel(np.array([imp_fit, imp_low, imp_high]))
        plt_dat = plt_dat.append(mod_dat)
    return plt_dat

def prep_dat_cat_warming(impacts, loc, return_period, wl=np.array([0.7, 1.2, 1.5, 2.0]),
                         warming=np.array([str(2000), 'Current', '1.5 degree', '2 degree']), RP=100.):
    plt_dat = pd.DataFrame(columns=['Model', 'Warming', 'RR', 'Value'])
    
    for w, wlevel in enumerate(wl):
        for model, m_imp in impacts.items():
            rp_mod = return_period[model]
            imp = m_imp[loc][str(wlevel)]
            
            imp_fit = np.interp(RP, rp_mod, imp[0,:])
            imp_low = np.interp(RP, rp_mod, imp[1,:])
            imp_high = np.interp(RP, rp_mod, imp[2,:])
                       
            mod_dat = pd.DataFrame(columns=['Model', 'Warming', 'RR', 'Value'])
            mod_dat['Model'] = np.repeat(model, 3)
            mod_dat['RR'] = np.array(['fit', 'low', 'high'])
            mod_dat['Value'] = np.ravel(np.array([imp_fit, imp_low, imp_high]))
            mod_dat['Warming'] = np.tile(warming[w], 3)
            plt_dat = plt_dat.append(mod_dat)
            
    return plt_dat

def prep_dat_cat_changes_warming(changes, loc, cities, wl=np.array([1.2, 1.5, 2.0]),
                         warming=np.array(['Current', '1.5 degree', '2 degree'])):
    plt_dat = pd.DataFrame(columns=['Model', 'Warming', 'Value'])
    
    for w, wlevel in enumerate(wl):
        for model, m_change in changes.items():                       
            mod_dat = pd.DataFrame(columns=['Model', 'Warming', 'Value'])
            mod_dat['Model'] = np.repeat(model, 1)
            mod_dat['Warming'] = np.repeat(warming[w], 1)
            mod_dat['Value'] = np.repeat(m_change[wlevel][cities.index(loc)], 1)
            plt_dat = plt_dat.append(mod_dat)
            
    return plt_dat

def plot_cat_RP(impacts, return_period, locations, location_names, wl=2.0):
    dat = []
    for loc in locations:
        dat.append(prep_dat_cat_RP(impacts, loc, return_period, wl=wl))
                   
    for i, d in enumerate(dat):
        sns.set_style('white')
        sns.set_context(rc={"font.size":30,"axes.titlesize":36,
                            "axes.labelsize":30, "xtick.labelsize":24,
                            "ytick.labelsize":24}) 
        g = sns.catplot(x="Return_period", y="Value", hue="Model",
                      data=d, dodge=0.5, kind="point", linestyles='',
                      err_style='bar', errwidth=4.0, markers="+", scale=2.5,
                      palette='viridis',aspect=1.8/1, ci=100,
                      legend=False)
        g.set(xlabel="Return period", ylabel="Heat mortality [#]", ylim=(0,None),
              title=location_names[i])
        g.savefig('g'+str(i)+'.png')
        plt.close(g.fig)
    
    fig, axes = plt.subplots(4, 3, figsize = (24, 20))
    k = 0
    for i in range(4):
        for j in range(3):
            axes[i,j].imshow(mpimg.imread('g'+str(k)+'.png'))    
            axes[i,j].axis('off')
            k=k+1
            
def plot_cat_warming(impacts, return_period, locations, location_names):
    dat = []
    for loc in locations:
        dat.append(prep_dat_cat_warming(impacts, loc, return_period))
                   
    for i, d in enumerate(dat):
        sns.set_style('white')
        sns.set_context(rc={"font.size":30,"axes.titlesize":36,
                            "axes.labelsize":30, "xtick.labelsize":24,
                            "ytick.labelsize":24}) 
        g = sns.catplot(x="Warming", y="Value", hue="Model",
                      data=d, dodge=0.5, kind="point", linestyles='',
                      err_style='bar', errwidth=4.0, markers="+", scale=2.5,
                      palette='viridis',aspect=1.8/1, ci=100,
                      legend=False)
        g.set(xlabel="Warming level", ylabel="Heat mortality [#]", ylim=(0,None),
              title=location_names[i])
        g.savefig('g'+str(i)+'.png')
        plt.close(g.fig)
    
    fig, axes = plt.subplots(4, 3, figsize = (24, 20))
    k = 0
    for i in range(4):
        for j in range(3):
            axes[i,j].imshow(mpimg.imread('g'+str(k)+'.png'))    
            axes[i,j].axis('off')
            k=k+1
            
def plot_cat_warming_changes(changes, locations, cities, location_names,
                             ylab='', prc=False):
    dat = []
    for loc in locations:
        dat.append(prep_dat_cat_changes_warming(changes, loc, cities))
                       
    for i, d in enumerate(dat):
        if prc:
            d.Value = d.Value*100
        sns.set_style('white')
        sns.set_context(rc={"font.size":30,"axes.titlesize":36,
                            "axes.labelsize":30, "xtick.labelsize":24,
                            "ytick.labelsize":24}) 
        g = sns.catplot(x="Warming", y="Value", hue="Model",
                      data=d, dodge=0.5, kind="point", linestyles='',
                      err_style='bar', errwidth=4.0, markers="+", scale=2.5,
                      palette='viridis',aspect=1.8/1, ci=100,
                      legend=False)
        g.set(xlabel="Warming level", ylabel=ylab, ylim=(0,None),
              title=location_names[i])
        g.savefig('g'+str(i)+'.png')
        plt.close(g.fig)
    
    fig, axes = plt.subplots(4, 3, figsize = (24, 20))
    k = 0
    for i in range(4):
        for j in range(3):
            axes[i,j].imshow(mpimg.imread('g'+str(k)+'.png'))    
            axes[i,j].axis('off')
            k=k+1
            

# aa = prep_dat_cat_changes_warming(RP_changes, loc, cities)

# plot_cat_RP(impacts, return_period, locations, loc_names, wl=2.0)
# plot_cat_warming(impacts, return_period, locations, loc_names)
# plot_cat_warming_changes(RP_changes, locations, cities, loc_names,
#                              ylab='New return period')
# plot_cat_warming_changes(Imp_changes_pct, locations, cities, loc_names,
#                              ylab='Increase in mortality [%-diff]')
# plot_cat_warming_changes(Imp_changes_pct, locations, cities, loc_names,
#                              ylab='Change in mortality [%]', prc=True)


'''------------------------------------------------------------------------'''
''' bootstrap heat mortality reference period '''

def plot_ref_bootstrap(city, f_path_ref, imp_obs, rp_obs, all_rr, all_mort,
                       all_tmm, n_sample=1000, conf_int=0.95, col='mediumorchid',
                       tickfont=20, labelfont=22, ylab=True):
    
    data_temp_ref = pd.read_hdf(f_path_ref, key=city)
    data_temp_ref = data_temp_ref.iloc[0:(data_temp_ref.shape[0]//365)*365,:]
    RR = pd.read_excel(all_rr, sheet_name=city)
    mort = pd.read_excel(all_mort, sheet_name=city)
    tmm = all_tmm.TMM[all_tmm.city==city].values[0]

    mort_imp, rp = calc_heat_mortality(data_temp_ref, mort, RR, tmm)
    impact, return_per = calc_bootstrapping_ImpFreqCurve(mort_imp[0], n_sample=data_temp_ref.shape[0]//365, conf_int=0.95)
    
    plt.plot(return_per, impact[:, 0], color=col, linewidth=3)
    plt.fill_between(return_per, impact[:, 1], impact[:, 2],
                        color=col, alpha=0.15)
    plt.plot(rp_obs[city], imp_obs[city][0], color='black', linewidth=3)
    plt.xscale('log')
    plt.xticks(fontsize=tickfont)
    plt.yticks(fontsize=tickfont)
    plt.xlim(1,len(rp_obs[city]))
    #plt.ylim(0,False)
    plt.xlabel("Return period (Year)", fontsize=labelfont)
    if ylab:
        plt.ylabel("Heat mortality (#)", fontsize=labelfont)
    sns.despine()
    
    return impact, return_per

def plot_grid_line_ref_bootstrap(dat_path, models, locations, imp_obs, rp_obs, all_rr, all_mort, all_tmm):
    
    plt.figure(figsize=(20, 24))
    #plt.suptitle(title, fontsize=44)
    cm = plt.get_cmap('twilight')
    cols = [cm(i / (len(models)+1)) for i in range(len(models)+1)]
    i = 0
    for m, model in enumerate(models):
        file = dat_path.joinpath('reference_data_BC_'+model+'.h5')
        for l, loc in enumerate(locations):
            plt.subplot(len(models),len(locations),i+1)
            plt.subplots_adjust(hspace=.4)
            plt.subplots_adjust(wspace=.4)
            if l==0:
                plot_ref_bootstrap(loc, file, imp_obs, rp_obs, all_rr, all_mort,
                                       all_tmm, col=cols[m+1])
                plt.text(0.25, 0.2, model, rotation=90, fontsize=24)
            else:
                plot_ref_bootstrap(loc, file, imp_obs, rp_obs, all_rr, all_mort,
                                       all_tmm, col=cols[m+1], ylab=False)
            if m==0:
                plt.title(loc)

            i = i+1
            

# dat_path= DP.joinpath('ESM/Global/run_220216/')
# sample_loc =['spal.bra9718','toky.jap7215','zrch.sui9513']
# plot_grid_line_ref_bootstrap(dat_path, LE_models, sample_loc, imp_obs, rp_obs, all_rr, all_mort, all_tmm)


'''------------------------------------------------------------------------'''
''' Plot petridish '''

class Point:
    """A little class representing an SVG circle."""

    def __init__(self, cx, cy, r, icolour=None):
        """Initialize the circle with its centre, (cx,cy) and radius, r.

        icolour is the index of the circle's colour.

        """
        self.cx, self.cy, self.r = cx, cy, r
        self.icolour = icolour

    def overlap_with(self, cx, cy, r):
        """Does the circle overlap with another of radius r at (cx, cy)?"""

        d = np.hypot(cx-self.cx, cy-self.cy)
        return d < r + self.r

class PointCloud:
    """A class for drawing circles-inside-a-circle.
    Inspired by https://scipython.com/blog/packing-circles-in-a-circle/
    """

    def __init__(self, width=600, height=600, R=250, n=800, rho_min=0.005,
                 rho_max=0.05, colours=None):
        """Initialize the Circles object.

        width, height are the SVG canvas dimensions
        R is the radius of the large circle within which the small circles are
        to fit.
        n is the maximum number of circles to pack inside the large circle.
        rho_min is rmin/R, giving the minimum packing circle radius.
        rho_max is rmax/R, giving the maximum packing circle radius.
        colours is a list of SVG fill colour specifiers to be referenced by
            the class identifiers c<i>. If None, a default palette is set.

        """

        self.width, self.height = width, height
        self.R, self.n = R, n
        # The centre of the canvas
        self.CX, self.CY = self.width // 2, self.height // 2
        self.rmin, self.rmax = R * rho_min, R * rho_max
        self.colours = colours or ['#993300', '#a5c916', '#00AA66', '#FF9900']

   
    def _place_circle(self, r):
        # The guard number: if we don't place a circle within this number
        # of trials, we give up.
        guard = 500
        while guard:
            # Pick a random position, uniformly on the larger circle's interior
            cr, cphi = ( self.R * np.sqrt(np.random.random()),
                         2*np.pi * np.random.random() )
            cx, cy = cr * np.cos(cphi), cr * np.sin(cphi)
            if cr+r < self.R:
            # The circle fits inside the larger circle.
                if not any(circle.overlap_with(self.CX+cx, self.CY+cy, r)
                                    for circle in self.circles):
                    # The circle doesn't overlap any other circle: place it.
                    circle = Point(cx+self.CX, cy+self.CY, r,
                                icolour=np.random.randint(len(self.colours)))
                    self.circles.append(circle)
                    return
            guard -= 1
        # Warn that we reached the guard number of attempts and gave up for
        # for this circle.
        print('guard reached.')

    def make_circles(self):
        """Place the little circles inside the big one."""

        # First choose a set of n random radii and sort them. We use
        # random.random() * random.random() to favour small circles.
        self.circles = []
        r = self.rmin + (self.rmax - self.rmin) * np.random.random(
                                self.n) * np.random.random(self.n)
        r[::-1].sort()
        # Do our best to place the circles, larger ones first.
        for i in range(self.n):
            self._place_circle(r[i])
            
def plot_petri_dish(impacts, return_period, location, model='Median', title='',
                    sort_points=False):
    
    wl=np.array([1.2, 1.5, 2.0])
    wl_name = ['Risk 2000', 'Risk 2020', 'Risk 1.5$^\circ$C', 'Risk 2.0$^\circ$C']  

    col_lvls, al_lvls = prep_data_petri_dish(impacts, return_period,
                                                 location, model, wl)
    
    fig, axes = plt.subplots(1, 4, figsize = (24, 6))
    plt.rcParams.update({'font.family':'arial'})
    plt.suptitle(title, fontsize=44)
    

    for w, wlevel in enumerate(wl_name):
        cloud = PointCloud(width=100, height=100, R=50, n=100, rho_min=0.065,
                     rho_max=0.065)
        cloud.make_circles()
        while len(cloud.circles)<100:
            cloud = PointCloud(width=100, height=100, R=50, n=100, rho_min=0.065,
                         rho_max=0.065)
            cloud.make_circles()
        cx = np.zeros(100)
        cy = np.zeros(100)
        for j in range(100):
            cx[j] = cloud.circles[j].cx
            cy[j] = cloud.circles[j].cy
        if sort_points:
            r = (cx-50.)**2+(cy-50)**2
            r_ind = 99-np.argsort(r.argsort())
            col = []
            al = []
            for i in r_ind:
                col.append(col_lvls[w][i])
                al.append(al_lvls[w][i])
        else:
            col = col_lvls[w]
            al = al_lvls[w]
        circle = Circle((50, 50), 50, facecolor='none',
                        edgecolor='black', linewidth=1, alpha=1)
    
        axes[w].scatter(cx, cy, s=400, marker='o', c=col, alpha=al)
        axes[w].add_patch(circle)
        axes[w].get_xaxis().set_visible(False)
        axes[w].get_yaxis().set_visible(False)
        axes[w].axis('off')
        axes[w].set_title(wlevel, fontsize=32)   
    plt.subplots_adjust(left=0.15,
                    bottom=0.1, 
                    right=0.85, 
                    top=0.8, 
                    wspace=0.1, 
                    hspace=0.1)
    
    cm = plt.get_cmap('inferno_r')
    cols = [cm(i / 10) for i in range(10)]
    
    # p1 = plt.scatter(np.NaN, np.NaN, color=cols[1], s=400, marker='o')
    p2 = plt.plot(np.NaN, np.NaN, color=cols[3], markersize=20, marker='o', linestyle='')
    p3 = plt.plot(np.NaN, np.NaN, color=cols[6], markersize=20, marker='o', linestyle='')
    p4 = plt.plot(np.NaN, np.NaN, color=cols[9], markersize=20, marker='o', linestyle='')



    plt.legend([p2[0],p3[0], p4[0]],
               ['1-in-10 year season in 2000',
                '1-in-100-year season in 2000',
                '1-in-500-year season in 2000'],
               bbox_to_anchor=(-1.1, -0.25), loc='lower center', ncol=3,
               fontsize=24, frameon=False) 
    plt.show()


def prep_data_petri_dish(impacts, return_period, location, model, wl):
    cm = plt.get_cmap('inferno_r')
    cols = [cm(i / 10) for i in range(10)]

    # make colors for petri dishes
    c = np.arange(0,100,1); col_lvls = dict(); al_lvls = dict()

    # base
    al = []; col = []
    for j in c:
        if j < 90: col.append(cols[1]); al.append(0.4)
        elif j < 99: col.append(cols[3]); al.append(0.6)
        else: col.append(cols[6]); al.append(0.8)
    col_lvls[0] = col
    al_lvls[0] = al

    for w, wlevel in enumerate(wl):
        c_10 = change_in_RP(10, return_period[model], impacts[model][location][str(0.7)][0,:],
                                              return_period[model], impacts[model][location][str(wlevel)][0,:])
        c_100 = change_in_RP(100, return_period[model], impacts[model][location][str(0.7)][0,:],
                                               return_period[model], impacts[model][location][str(wlevel)][0,:])
        c_500 = change_in_RP(500, return_period[model],impacts[model][location][str(0.7)][0,:],
                                                return_period[model], impacts[model][location][str(wlevel)][0,:])
        al = []; col = []
        for j in c:
            if j < 100-np.round(1/c_10*100): col.append(cols[1]); al.append(0.4)
            elif j < 100-np.round(1/c_100*100): col.append(cols[3]); al.append(0.6)
            elif j < 100-np.round(1/c_500*100): col.append(cols[6]); al.append(0.8)
            else: col.append(cols[9]); al.append(1)
    
        col_lvls[w+1] = col
        al_lvls[w+1] = al
        
    return col_lvls, al_lvls



# plot_petri_dish(impacts, return_period, location=locations[10], model='Median', title='')
# plot_petri_dish(impacts, return_period, location=locations[10], model='Median', title='', sort_points=True)

# for c in locations:
#     plot_petri_dish(impacts, return_period, location=c, model='Median', title=c, sort_points=True)

''' Appendix model change '''

def plot_grid_figure_model_change(locations, loc_names, RP_changes, pct_changes,
                                  ylim_RP=np.nan, ylim_pct=np.nan):

    plt.figure(figsize=(16, 10))
    grid = plt.GridSpec(2, len(locations), wspace=0.3, hspace=0.6) # create a grid for the subplots
    
    for l, loc in enumerate(locations):
        
        if l==0: ylab=True
        else: ylab=False

        if not any(np.isnan(ylim_RP)): ymax=ylim_RP[l]
        else: ymax = 0
        ax = plt.subplot(grid[0,l])
        plot_model_RP_changes(loc, RP_changes, title=loc_names[l],
                            cmap='twilight', ylab=ylab, ymax=ymax, labelfont=24)
        
        if not any(np.isnan(ylim_pct)): ymax=ylim_pct[l]
        else: ymax = 0
        ax = plt.subplot(grid[1,l])
        plot_model_pct_changes(loc, pct_changes, title='',
                            cmap='twilight', ylab=ylab, ymax=ymax, labelfont=24)
    
        # legends
        if l==1:
            cm = plt.get_cmap('twilight')
            cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
            p1 = plt.scatter(np.NaN, np.NaN, marker='d', s=100, color=cols[1])
            p2 = plt.scatter(np.NaN, np.NaN, marker='d', s=100, color=cols[2])
            p3 = plt.scatter(np.NaN, np.NaN, marker='d', s=100, color=cols[3])
            p4 = plt.scatter(np.NaN, np.NaN, marker='d', s=100, color=cols[4])
            p5 = plt.scatter(np.NaN, np.NaN, marker='d', s=100, color=cols[5])
    
            plt.legend([p1, p2, p3, p4, p5],
                       LE_models,
                       bbox_to_anchor=(0.5, -0.4), loc='lower center', ncol=5,
                       fontsize=18, frameon=False)

# plot_grid_figure_model_change(locloc, locloc, RP_changes, Pct_mort,
#                               ylim_RP=[60,50,25], ylim_pct=[6,24,18])


def plot_model_RP_changes(loc, changes, cities=cities, wl=[1.2, 1.5, 2.0],
                       x_ticks=['2020', '1.5$^\circ$C', '2.0$^\circ$C'],
                       title='', cmap='viridis', LE_models=LE_models,
                       labelfont=22, titlefont=22, ylab=True, ymax=0,
                       y_label='New return period (y)', x_label=''):
    
    cm = plt.get_cmap(cmap)
    cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
    
    i_loc = cities.index(loc)
    x = np.arange(0, len(wl))+0.5
    for m, model in enumerate(LE_models):
        y = np.zeros([len(wl), 1])
        for w, warm in enumerate(wl):
            y[w] = changes[model][warm][i_loc]
        plt.scatter(x, y, marker='d',color=cols[m+1],s=100)
        
    
    if ylab:
        plt.ylabel(y_label, fontsize=labelfont)
    if ymax>0:
        plt.ylim(0,ymax)
    plt.xticks(x, x_ticks, fontsize=labelfont)
    plt.xlabel(x_label, fontsize=labelfont)
    plt.xlim(0,np.max(x)+0.5)
    plt.yticks(fontsize=labelfont)
    sns.despine()
    plt.title(title, fontsize=titlefont)
    
def plot_model_pct_changes(loc, changes, cities=cities, wl=[1.2, 1.5, 2.0],
                       x_ticks=['2020', '1.5$^\circ$C', '2.0$^\circ$C'],
                       title='', cmap='viridis', LE_models=LE_models,
                       labelfont=22, titlefont=22, ylab=True, ymax=0,
                       y_label='Heat mortality (%)', x_label=''):
    
    cm = plt.get_cmap(cmap)
    cols = [cm(i / (len(LE_models)+1)) for i in range(len(LE_models)+1)]
    
    i_loc = cities.index(loc)
    x = np.arange(0, len(wl))+0.5
    for m, model in enumerate(LE_models):
        y = np.zeros([len(wl), 1])
        for w, warm in enumerate(wl):
            y[w] = changes[model][warm][i_loc]*100
        plt.scatter(x, y, marker='d',color=cols[m+1],s=100)
        
    
    if ylab:
        plt.ylabel(y_label, fontsize=labelfont)
    if ymax>0:
        plt.ylim(0,ymax)
    #plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter())
    plt.xticks(x, x_ticks, fontsize=labelfont)
    plt.xlabel(x_label, fontsize=labelfont)
    plt.xlim(0,np.max(x)+0.5)
    plt.yticks(fontsize=labelfont)
    sns.despine()
    plt.title(title, fontsize=titlefont)

