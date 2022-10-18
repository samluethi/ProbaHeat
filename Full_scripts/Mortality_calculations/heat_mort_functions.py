"""
Created on Oct 17 2022
@author: Samuel Luethi - samuel.luethi@usys.ethz.ch

"""

import numpy as np
import pandas as pd
import datetime as dt
from sklearn.utils import resample

# This script contains all functions for the heat mortality calculations.
# See script Calc_heat_mort.py for running the functions.


def calc_attrib_deaths_ens(ens, mort_doy, IF_RR, IF_temp, TMM):
    
    """ Calculates attributubal deaths to temperature and heat.

    Parameters
    ----------
    ens : np.array
        array containing temperature data
    mort_doy : np.array
        daily mortality data [n_years, 365]
    IF_RR : np.array
        relative risk of mortality
    IF_temp : np.array
        temperature axis corresponding to IF_RR
    TMM : float
        temperature of minimum mortality

    Returns
    -------
    atrib_deaths_temp : np.array
        daily temperature related deaths
    atrib_deaths_heat : np.array
        daily heat related deaths

    """
    
    n_years = int((ens.shape[0])/365)
    try: n_ens = ens.shape[1]
    except IndexError: n_ens = 1
    
    atrib_frac_temp = np.interp(ens, IF_temp, IF_RR)-1
    if n_ens > 1:
        atrib_deaths_temp = atrib_frac_temp * \
            np.transpose(np.tile(mort_doy, (n_ens, n_years)))
    else:
        atrib_deaths_temp = atrib_frac_temp * \
            np.tile(mort_doy, n_years)
            
    atrib_deaths_heat = atrib_deaths_temp.copy()
    atrib_deaths_heat[TMM>ens] = 0

    return atrib_deaths_temp, atrib_deaths_heat

def calc_imp_year(ens_mort):

    """ Summarizes daily mortality data to anual.

    Parameters
    ----------
    ens_mort : np.array
        array containing daily mortality

    Returns
    -------
    annual_mort : np.array
        annual deaths

    """

    n_years = int((ens_mort.shape[0])/365)
    try: n_ens = ens_mort.shape[1]
    except IndexError: n_ens = 1

    annual_mort = np.zeros([n_years, n_ens])
    if n_ens > 1:
        for i in range(n_years):
            annual_mort[i,:] = np.nansum(ens_mort[i*365:(i+1)*365,:], axis=0)
    else:
        for i in range(n_years):
            annual_mort[i] = np.nansum(ens_mort[i*365:(i+1)*365], axis=0)

    return annual_mort

def calc_ImpFreqCurve(mort_imp, freq=None):

    """ Calclate impact frequency

    Parameters
    ----------
    mort_imp : np.array
        array containing annual mortality

    Returns
    -------
    impact : np.array
        annual deaths
    return_per : np.array
        frequncy corresponding to impact

    """

    event = mort_imp.ravel()
    if freq==None:
        freq = np.zeros(len(event))
        freq[:] = 1./(len(event))
    sort_idxs = np.argsort(event)[::-1]
    # Calculate exceedence frequency
    exceed_freq = np.cumsum(freq[sort_idxs])
    # Set return period and imact exceeding frequency
    return_per = 1 / exceed_freq[::-1]
    impact = event[sort_idxs][::-1]
    
    return impact, return_per

def calc_bootstrapping_ImpFreqCurve(mort_imp, n_sample=0, sample_size=1000, conf_int=0.95):

    """ Calclates uncertainty of impact frequency using bootstrapping.

    Parameters
    ----------
    mort_imp : np.array
        array containing annual mortality
    n_sample : int
        number of samples of bootstrapping
    conf_int : float
        confidence interval

    Returns
    -------
    impact : np.array
        annual deaths [modelled, lower_bound, upper_bound]
    return_per : np.array
        frequncy corresponding to impact

    """
    # array to store impact frequency curves
    I = np.zeros([n_sample, sample_size])
    
    # bootstrapping
    if n_sample==0:
        n_sample = len(mort_imp.ravel())
    for i in range(sample_size):
        x = resample(mort_imp.ravel(), replace=True, n_samples=n_sample)
        ifc, rp = calc_ImpFreqCurve(x)
        I[:,i] = ifc
        
    # get confidence interval
    impact = np.zeros([n_sample, 3])   
    conf = (1-conf_int)/2
    impact[:,0] = np.percentile(I, 50, axis=1)
    impact[:,1] = np.percentile(I, conf*100, axis=1)
    impact[:,2] = np.percentile(I, (1-conf)*100, axis=1)
    
    return impact, rp

def wrapper_heat_mortality(ens, ens_ref, mort, RR):
    
    # process temperarure data
    # ens_ref = _remove_leap_days(ens_ref, ens_ref.time)
    ens_ref_T_only = ens_ref.drop(['lat', 'lon', 'time'], axis=1)
    # ens = _remove_leap_days(ens, ens.time)
    ens_T_only = ens.drop(['lat', 'lon', 'time'], axis=1)
    T = bias_correct_ensemble(ens_ref_T_only, mort.temp, ens_T_only)

    # define temperature of minimum mortality
    TMM = RR.temp.values[RR.RRfit.argmin()]
    
    # process mortality data
    # make annual series of mean daily mortality
    # remove leap days
    mort_clean = _remove_leap_days(mort, mort.date)
    # make annual series of mean daily mortality
    n_years = int(len(mort_clean.deaths.values)/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_doy_mean = np.mean(mort_doy,0)
    
    # calc attributable deaths
    ad_temp, ad_heat = calc_attrib_deaths_ens(T, mort_doy_mean,
                                        RR.RRfit.values, RR.temp.values, TMM)
    ad_temp_high, ad_heat_high = calc_attrib_deaths_ens(T, mort_doy_mean,
                                        RR.RRhigh.values, RR.temp.values, TMM)
    ad_temp_low, ad_heat_low = calc_attrib_deaths_ens(T, mort_doy_mean,
                                        RR.RRlow.values, RR.temp.values, TMM)
    # summarize to years
    annual_mort = calc_imp_year(ad_heat)
    annual_mort_high = calc_imp_year(ad_heat_high)
    annual_mort_low = calc_imp_year(ad_heat_low)
    
    # calc impact frequncy incl. uncertainty
    impact, return_period = calc_bootstrapping_ImpFreqCurve(annual_mort)
    impact_high, _ = calc_bootstrapping_ImpFreqCurve(annual_mort_high)
    impact_low, _ = calc_bootstrapping_ImpFreqCurve(annual_mort_low)
    
    return [impact, impact_low, impact_high], return_period

def wrapper_heat_mortality_obs(mort, RR):
    
    # define temperature of minimum mortality
    TMM = RR.temp.values[RR.RRfit.argmin()]
    
    # process mortality data
    # make annual series of mean daily mortality
    # remove leap days
    mort_clean = _remove_leap_days(mort, mort.date)
    # make annual series of mean daily mortality
    n_years = int(len(mort_clean.deaths.values)/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_doy_mean = np.mean(mort_doy,0)
    
    # calc attributable deaths
    ad_temp, ad_heat = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRfit.values, RR.temp.values, TMM)
    ad_temp_high, ad_heat_high = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRhigh.values, RR.temp.values, TMM)
    ad_temp_low, ad_heat_low = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRlow.values, RR.temp.values, TMM)
    # summarize to years
    annual_mort = calc_imp_year(ad_heat)
    annual_mort_high = calc_imp_year(ad_heat_high)
    annual_mort_low = calc_imp_year(ad_heat_low)
    
    # calc impact frequncy incl. uncertainty
    impact, return_period = calc_ImpFreqCurve(annual_mort)
    impact_high, _ = calc_ImpFreqCurve(annual_mort_high)
    impact_low, _ = calc_ImpFreqCurve(annual_mort_low)
    
    return [impact, impact_low, impact_high], return_period

def _remove_leap_days(data):
    date_series = pd.to_datetime(data.date)
    idx_leap_day = np.where((date_series.dt.is_leap_year.values==True) & \
                            (date_series.dt.day_of_year.values==60))[0]
    data_clean = data.drop(idx_leap_day)
    return data_clean

def _clean_data_to_365day_year(data):
    ''' This function adds nan values to a time series data frame in order to
    fill up values to full years and removes leap days '''
    
    d = pd.DataFrame()
    d['date'] = pd.date_range(start=dt.datetime(data.date.min().year, 1 , 1),
                             end=dt.datetime(data.date.max().year, 12 , 31), freq='D')
        
    d = d.merge(data, on='date', how='left')
    # correct for leap days
    dat = _remove_leap_days(d)
    
    return dat

def change_in_RP(RP, r_prior, i_prior, r_post, i_post):
    imp_level = np.interp(RP, r_prior, i_prior)
    imp_RP_post = np.interp(imp_level, i_post, r_post)
    
    return imp_RP_post

def change_in_mort(RP, impact_prior, impact_post, return_period):
    i_pr = np.interp(RP, return_period, impact_prior)
    i_po = np.interp(RP, return_period, impact_post)
    rel_dif = i_po/i_pr-1
    
    return rel_dif

def change_in_mort_pct(RP, impact_prior, impact_post, return_period, mort):

    mort_clean = _clean_data_to_365day_year(mort)
    # make annual series of mean daily mortality
    n_years = int(mort_clean.shape[0]/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_doy_mean = np.nanmean(mort_doy,0)
    annual_mort = np.nansum(mort_doy_mean)
    
    i_pr = np.interp(RP, return_period, impact_prior)
    i_po = np.interp(RP, return_period, impact_post)
    
    rel_heat_mort_prio = i_pr/annual_mort
    rel_heat_mort_post = i_po/annual_mort

    
    return rel_heat_mort_post-rel_heat_mort_prio

def pct_heat_mort(RP, impact, return_period, mort):
    
    mort_clean = _clean_data_to_365day_year(mort)
    # make annual series of mean daily mortality
    n_years = int(mort_clean.shape[0]/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_doy_mean = np.nanmean(mort_doy,0)
    annual_mort = np.nansum(mort_doy_mean)
    
    imp = np.interp(RP, return_period, impact)
    rel_heat_mort = imp/annual_mort
    
    
    return rel_heat_mort

def calc_heat_mortality_obs(mort, RR, TMM):
    
    """ Calculates annual heat mortality.

    Parameters
    ----------
    T : pd.DataFrame
        array containing bias corrected temperature data
    mort : pd.DataFrame
        df containing date, daily mortality and temperature data
    RR : pd.DataFrame
        relative risk of mortality given temperature

    Returns
    -------
    impact : np.array
        annual deaths [modelled, lower_bound, upper_bound]
    return_per : np.array
        frequncy corresponding to impact

    """
    
    # clean mortality data
    mort_clean = _clean_data_to_365day_year(mort)
    
    # make annual series of mean daily mortality
    n_years = int(mort_clean.shape[0]/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_doy_mean = np.nanmean(mort_doy,0)
    
    # calc attributable deaths
    ad_temp, ad_heat = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRfit.values, RR.temp.values, TMM)
    ad_temp_high, ad_heat_high = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRhigh.values, RR.temp.values, TMM)
    ad_temp_low, ad_heat_low = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRlow.values, RR.temp.values, TMM)

    # summarize to years
    annual_mort = calc_imp_year(ad_heat)
    annual_mort_high = calc_imp_year(ad_heat_high)
    annual_mort_low = calc_imp_year(ad_heat_low)
    
    # calc impact frequncy incl. uncertainty
    impact, return_period = calc_ImpFreqCurve(annual_mort)
    impact_high, _ = calc_ImpFreqCurve(annual_mort_high)
    impact_low, _ = calc_ImpFreqCurve(annual_mort_low)
    
    return np.vstack((impact, impact_low, impact_high)), return_period

def calc_heat_mortality_fraction_obs(mort, RR, TMM):
    
    """ Calculates annual heat mortality.

    Parameters
    ----------
    T : pd.DataFrame
        array containing bias corrected temperature data
    mort : pd.DataFrame
        df containing date, daily mortality and temperature data
    RR : pd.DataFrame
        relative risk of mortality given temperature

    Returns
    -------
    impact : np.array
        annual deaths [modelled, lower_bound, upper_bound]
    return_per : np.array
        frequncy corresponding to impact

    """
    
    # clean mortality data
    mort_clean = _clean_data_to_365day_year(mort)
    
    # make annual series of mean daily mortality
    n_years = int(mort_clean.shape[0]/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_year = np.nansum(mort_doy, 1)
    mort_doy_mean = np.nanmean(mort_doy,0)
    
    # calc attributable deaths
    ad_temp, ad_heat = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRfit.values, RR.temp.values, TMM)
    ad_temp_high, ad_heat_high = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRhigh.values, RR.temp.values, TMM)
    ad_temp_low, ad_heat_low = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRlow.values, RR.temp.values, TMM)

    # summarize to years
    annual_mort = calc_imp_year(ad_heat)
    annual_mort_high = calc_imp_year(ad_heat_high)
    annual_mort_low = calc_imp_year(ad_heat_low)
    
    a = np.zeros(len(annual_mort))
    a_high = np.zeros(len(annual_mort))
    a_low = np.zeros(len(annual_mort))
    for i in range(len(a)):
        a[i] = annual_mort[i][0]
        a_high[i] = annual_mort_high[i][0]
        a_low[i] = annual_mort_low[i][0]
    # calc impact frequncy incl. uncertainty
    impact, return_period = calc_ImpFreqCurve(a/mort_year)
    impact_high, _ = calc_ImpFreqCurve(a_high/mort_year)
    impact_low, _ = calc_ImpFreqCurve(a_low/mort_year)
    
    return np.vstack((impact, impact_low, impact_high)), return_period

def calc_heat_mortality(T, mort, RR, TMM):
    
    """ Calculates annual heat mortality.

    Parameters
    ----------
    T : pd.DataFrame
        array containing bias corrected temperature data
    mort : pd.DataFrame
        df containing date, daily mortality and temperature data
    RR : pd.DataFrame
        relative risk of mortality given temperature

    Returns
    -------
    impact : np.array
        annual deaths [modelled, lower_bound, upper_bound]
    return_per : np.array
        frequncy corresponding to impact

    """
   
    # process mortality data
    # make annual series of mean daily mortality
    # clean mortality data
    mort_clean = _clean_data_to_365day_year(mort)
    
    # make annual series of mean daily mortality
    n_years = int(mort_clean.shape[0]/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_doy_mean = np.nanmean(mort_doy,0)
    
    # calc attributable deaths  
    ad_temp, ad_heat = calc_attrib_deaths_ens(T, mort_doy_mean,
                                        RR.RRfit.values, RR.temp.values, TMM)
    ad_temp_high, ad_heat_high = calc_attrib_deaths_ens(T, mort_doy_mean,
                                        RR.RRhigh.values, RR.temp.values, TMM)
    ad_temp_low, ad_heat_low = calc_attrib_deaths_ens(T, mort_doy_mean,
                                        RR.RRlow.values, RR.temp.values, TMM)
    # summarize to years
    annual_mort = calc_imp_year(ad_heat)
    annual_mort_high = calc_imp_year(ad_heat_high)
    annual_mort_low = calc_imp_year(ad_heat_low)
    
    # calc impact frequncy incl. uncertainty
    impact, return_period = calc_ImpFreqCurve(annual_mort)
    impact_high, _ = calc_ImpFreqCurve(annual_mort_high)
    impact_low, _ = calc_ImpFreqCurve(annual_mort_low)
    
    return np.vstack((impact, impact_low, impact_high)), return_period


def get_rank_heat_mort_obs(mort, RR, TMM):
    
    """ Calculates annual heat mortality.

    Parameters
    ----------
    T : pd.DataFrame
        array containing bias corrected temperature data
    mort : pd.DataFrame
        df containing date, daily mortality and temperature data
    RR : pd.DataFrame
        relative risk of mortality given temperature

    Returns
    -------
    impact : np.array
        annual deaths [modelled, lower_bound, upper_bound]
    return_per : np.array
        frequncy corresponding to impact

    """
    
    # clean mortality data
    mort_clean = _clean_data_to_365day_year(mort)
    
    # make annual series of mean daily mortality
    n_years = int(mort_clean.shape[0]/365)
    mort_doy = np.reshape(mort_clean.deaths.values, (n_years, 365))
    mort_doy_mean = np.nanmean(mort_doy,0)
    
    # calc attributable deaths
    ad_temp, ad_heat = calc_attrib_deaths_ens(mort_clean.temp, mort_doy_mean,
                                        RR.RRfit.values, RR.temp.values, TMM)
    
    # summarize to years
    annual_mort = calc_imp_year(ad_heat)
    a = np.zeros(len(annual_mort))
    for i in range(len(a)):
        a[i] = annual_mort[i][0]
    
    years = np.arange(mort.date.min().year, mort.date.max().year+1, 1)
    
    rank_year = years[np.argsort(a)]
    rank_imp = a[np.argsort(a)]
    
    return rank_imp, rank_year

def calc_haz_extremes(T, TMM):
    
    """ Calculates annual heat mortality.

    Parameters
    ----------
    T : pd.DataFrame
        array containing bias corrected temperature data
    TMM : np.float
        minimum mortality temperature
    Returns
    -------
    impact : np.array
        annual hazard
    return_per : np.array
        frequncy corresponding to impact

    """
    Tdat = T.to_numpy()
    Tdat = Tdat-TMM
    Tdat[0>Tdat] = 0
    # summarize to years
    annual_haz = calc_imp_year(Tdat)
    
    # calc impact frequncy incl. uncertainty
    impact, return_period = calc_ImpFreqCurve(annual_haz)
    
    return impact, return_period