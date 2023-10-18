'''
Module docstring
'''
import matplotlib
matplotlib.use('Agg')
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import os,sys
import pickle as pickle #cPickle is an optimised version of Pickle and is O(100) times faster
import glob
from collections import OrderedDict
import subprocess
import scipy
from scipy.stats.mstats import pearsonr

sys.path.append('/data/users/hadom/branches/git/eerie-project/storm_track_analysis/assess')
import tc_assess_code
import tc_assess_code.calc_interann
import ts_obs

BASIN_NAME = tc_assess_code.BASIN_NAME

import time

def lineplot_setup(xaxis_range,yaxis_range,title,xtitle,ytitle):
    """ 
    Set up a lineplot 
    
    """
#    fig.subplots_adjust(right=0.78, left=0.13, bottom=0.1, hspace=0.25)
    plt.grid()
    plt.title('%s' % title)  
    ax = plt.gca()
    ax.set_xlabel(xtitle, fontsize='large')
    ax.set_ylabel(ytitle, ha='center', color='black', fontsize='large')
    ax.set_xlim(xaxis_range[0]-1, xaxis_range[len(xaxis_range)-1]+1)
    ax.set_ylim(yaxis_range[0],yaxis_range[1])

def calculate_hurdat_frequency(obs_storms, annual_count_obs, obs_ds, BASINS, months_nh, months_sh):
    for ib, basin in enumerate(BASINS):      
        count = None
        count_cat = None
        if (BASINS[basin][0] == 'nh'):
            months = months_nh
        elif (BASINS[basin][0] == 'sh'):
            months = months_sh

        # calculate counts
        count = np.ma.asarray(tc_assess_code.annual_vmax_storm_counts(obs_storms, annual_count_obs['years'], months, basin, storm_types=['TS', 'HU', 'MH']), np.int64)
        annual_count_obs[obs_ds].update({('tc',basin): count})

        # calculate ace
        ace = np.ma.asarray(tc_assess_code.annual_vmax_storm_ace(obs_storms, annual_count_obs['years'], months, basin, storm_types=['TS', 'HU', 'MH']), np.int64)
        annual_count_obs[obs_ds].update({('tc', 'ace', basin): ace})

        # calculate counts in categories
        count_cat = np.ma.asarray(tc_assess_code.annual_vmax_storm_counts_category(obs_storms, annual_count_obs['years'], months, basin, storm_types=['TS', 'HU'
, 'MH']), np.int64)
        for icat in range(6):
            annual_count_obs[obs_ds].update({('tc',basin, icat): count_cat[:,icat]})

        # calculate count for hurricanes only
        count = np.ma.asarray(tc_assess_code.annual_vmax_storm_counts(obs_storms, annual_count_obs['years'], months, basin, storm_types=['HU', 'MH']), np.int64)
        annual_count_obs[obs_ds].update({('hur', basin) : count})

        # calculate ace for hurricanes only
        ace = np.ma.asarray(tc_assess_code.annual_vmax_storm_ace(obs_storms, annual_count_obs['years'], months, basin, storm_types=['HU', 'MH']), np.int64)
        annual_count_obs[obs_ds].update({('hur', 'ace', basin) : ace})

        count_cat = np.ma.asarray(tc_assess_code.annual_vmax_storm_counts_category(obs_storms, annual_count_obs['years'], months, basin, storm_types=['HU', 'MH'
]), np.int64)
        for icat in range(6):
            annual_count_obs[obs_ds].update({('hur', basin, icat): count_cat[:,icat]})

        ace_cat = np.ma.asarray(tc_assess_code.annual_vmax_storm_ace_category(obs_storms, annual_count_obs['years'], months, basin, storm_types=['TS', 'HU', 'MH'
]), np.int64)
        for icat in range(6):
            annual_count_obs[obs_ds].update({('tc', 'ace', basin, icat): ace_cat[:,icat]})

        ace_cat = np.ma.asarray(tc_assess_code.annual_vmax_storm_ace_category(obs_storms, annual_count_obs['years'], months, basin, storm_types=['HU', 'MH'
]), np.int64)
        for icat in range(6):
            annual_count_obs[obs_ds].update({('hur', 'ace', basin, icat): ace_cat[:,icat]})


    # data we use for some basins only extends to ~1970s, so need to mask it
    mask_obs_data(annual_count_obs, obs_ds)

def mask_obs_data(annual_count_obs, obs_ds):
    '''
    Mask out years where there is no good data in some basins
    '''
    # SH (including S Pacific and S Indian) - no good data before 1985
# now do SH
    if obs_ds =='hurdat2':
        miss = np.where(annual_count_obs['years'] <= 1984)
        for bas in ['sp', 'si', 'au', 'sh']:
    #print 'miss ',miss,type(miss)
            annual_count_obs[obs_ds][('tc', bas)][miss] = np.ma.masked
            annual_count_obs[obs_ds][('hur', bas)][miss] = np.ma.masked
            for icat in range(6):
                annual_count_obs[obs_ds][('hur',bas, icat)][miss] = np.ma.masked
                annual_count_obs[obs_ds][('tc',bas, icat)][miss] = np.ma.masked
    
        miss = np.where(annual_count_obs['years'] <= 1970)
        for bas in ['ni', 'nh']:
    #print 'miss ',miss,type(miss)
            annual_count_obs[obs_ds][('tc', bas)][miss] = np.ma.masked
            annual_count_obs[obs_ds][('hur', bas)][miss] = np.ma.masked
            for icat in range(6):
                annual_count_obs[obs_ds][('hur',bas, icat)][miss] = np.ma.masked
                annual_count_obs[obs_ds][('tc',bas, icat)][miss] = np.ma.masked

    elif obs_ds == 'ibtracs':
        miss = np.where(annual_count_obs['years'] <= 1981)
        for bas in BASINS_LIST:
            annual_count_obs[obs_ds][('tc', bas)][miss] = np.ma.masked
            annual_count_obs[obs_ds][('hur', bas)][miss] = np.ma.masked
            for icat in range(6):
                annual_count_obs[obs_ds][('hur',bas, icat)][miss] = np.ma.masked
                annual_count_obs[obs_ds][('tc',bas, icat)][miss] = np.ma.masked
            
def calculate_model_frequency(runid, runkey, annual_count_models, nensemble, storms_model, BASINS, months_nh, months_sh):
    
    for ib, basin in enumerate(BASINS):
        count = None
        count_category = None
        if (BASINS[basin][0] == 'nh'):
            months = months_nh
            years_modelf = annual_count_models[runkey][('years','nh')]
        elif (BASINS[basin][0] == 'sh'):
            months = months_sh
            years_modelf = annual_count_models[runkey][('years','sh')]

        print('process basin ',ib, basin)
        count = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_count(storms_model, years_modelf, months, basin, nensemble), np.float64)
        #print 'basin, count ',ib, basin, count
        if len(count) != len(years_modelf):
            print('warning, length of tc freq not equal to no years ',len(count), len(years_modelf))
        annual_count_models[runkey].update({basin: count})

        count_category = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_count_category(storms_model, years_modelf, months, basin, nensemble), np.float64)
        for icat in range(6):
            annual_count_models[runkey].update({(basin,icat): count_category[:,icat]})
        
# when the NH (or SH) total is zero, then this is missing, so set the appropriate hemisphere BASINS to missing
    miss = np.where(annual_count_models[runkey]['nh'] == 0)
    #print 'miss ',miss,type(miss)
    annual_count_models[runkey]['nh'][miss] = np.ma.masked
    model_years_nonzero = np.delete(annual_count_models[runkey][('years','nh')], miss[0])
    
    #print 'model years non-zero ',runkey, model_years_nonzero
    annual_count_models[runkey]['years','nh', 'nonzero']  = model_years_nonzero

    for icat in range(6):
        annual_count_models[runkey][('nh',icat)][miss] = np.ma.masked
    for ib, basin in enumerate(BASINS):
        if BASINS[basin][0] == 'nh':
            annual_count_models[runkey][basin][miss] = np.ma.masked
            for icat in range(6):
                annual_count_models[runkey][(basin,icat)][miss] = np.ma.masked
# now do SH
    miss = np.where(annual_count_models[runkey]['sh'] == 0)
    #print 'miss ',miss,type(miss)
    annual_count_models[runkey]['sh'][miss] = np.ma.masked
    model_years_nonzero = np.delete(annual_count_models[runkey][('years','sh')], miss[0])
    #print 'model SH years non-zero ',runkey, model_years_nonzero
    annual_count_models[runkey]['years','sh', 'nonzero'] = model_years_nonzero

    for icat in range(6):
        annual_count_models[runkey][('sh',icat)][miss] = np.ma.masked
    for ib, basin in enumerate(BASINS):
        if BASINS[basin][0] == 'sh':
            annual_count_models[runkey][basin][miss] = np.ma.masked
            for icat in range(6):
                annual_count_models[runkey][(basin,icat)][miss] = np.ma.masked
                
    return annual_count_models

def calculate_model_ace(runid, runkey, annual_count_models, nensemble, storms_model, BASINS, months_nh, months_sh):
    
    for ib, basin in enumerate(BASINS):
        ace = None
        if (BASINS[basin][0] == 'nh'):
            months = months_nh
            years_modelf = annual_count_models[runkey][('years','nh')]
        elif (BASINS[basin][0] == 'sh'):
            months = months_sh
            years_modelf = annual_count_models[runkey][('years','sh')]

        print('process basin ',ib, basin)
        ace = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_ace(storms_model, years_modelf, months, basin, nensemble), np.float64)
        #print('basin, ace ',ib, basin, ace)
        if len(ace) != len(years_modelf):
            print('warning, length of tc ace not equal to no years ',len(ace), len(years_modelf))

        annual_count_models[runkey].update({('ace', basin): ace})

        ace_category = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_ace_category(storms_model, years_modelf, months, basin, nensemble), np.float64)
        for icat in range(6):
            annual_count_models[runkey].update({('ace', basin, icat): ace_category[:,icat]})
        
# when the NH (or SH) total is zero, then this is missing, so set the appropriate hemisphere BASINS to missing
    miss = np.where(annual_count_models[runkey][('ace', 'nh')] == 0)
    #print 'miss ',miss,type(miss)
    annual_count_models[runkey][('ace', 'nh')][miss] = np.ma.masked

    for icat in range(6):
        annual_count_models[runkey][('ace', 'nh' ,icat)][miss] = np.ma.masked
    for ib, basin in enumerate(BASINS):
        if BASINS[basin][0] == 'nh':
            annual_count_models[runkey]['ace', basin][miss] = np.ma.masked
            for icat in range(6):
                annual_count_models[runkey][('ace', basin, icat)][miss] = np.ma.masked

# now do SH
    miss = np.where(annual_count_models[runkey][('ace', 'sh')] == 0)
    #print 'miss ',miss,type(miss)
    annual_count_models[runkey][('ace', 'sh')][miss] = np.ma.masked
    #model_years_nonzero = np.delete(annual_count_models[runkey][('years','sh')], miss[0])
    #print 'model SH years non-zero ',runkey, model_years_nonzero
    #annual_count_models[runkey]['years','sh', 'nonzero'] = model_years_nonzero

    for icat in range(6):
        annual_count_models[runkey][('ace','sh',icat)][miss] = np.ma.masked
    for ib, basin in enumerate(BASINS):
        if BASINS[basin][0] == 'sh':
            annual_count_models[runkey]['ace', basin][miss] = np.ma.masked
            for icat in range(6):
                annual_count_models[runkey][('ace',basin,icat)][miss] = np.ma.masked
                
    return annual_count_models

def calc_tc_frequency_mean_std(runs, BASINS, annual_count_obs, OBS_DATASETS, annual_count_models, OBS_TYPES, mean_count, std_count, \
                               mean_count_cat, std_count_cat):

    for ods in OBS_DATASETS:
        for ib, basin in enumerate(BASINS):
            for otype in OBS_TYPES:
                mean_count[(ods, otype, basin)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], basin)].mean(), dtype = np.float32)
                std_count[(ods, otype, basin)] =  np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype],basin)].std(), dtype = np.float32)
                for icat in np.arange(6):        
                    mean_count_cat[(ods, otype, basin, icat) ] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], basin, icat)].mean(), dtype = np.float32)
                    std_count_cat[(ods, otype, basin, icat)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], basin, icat)].std(), dtype = np.float32)

# do ACE
    for ods in OBS_DATASETS:
        for ib, basin in enumerate(BASINS):
            for otype in OBS_TYPES:
                mean_count[(ods, otype, 'ace', basin)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin)].mean(), dtype = np.float32)
                std_count[(ods, otype, 'ace', basin)] =  np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin)].std(), dtype = np.float32)
                for icat in np.arange(6):        
                    mean_count_cat[(ods, otype, 'ace', basin, icat) ] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin, icat)].mean(), dtype = np.float32)
                    std_count_cat[(ods, otype, 'ace', basin, icat)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin, icat)].std(), dtype = np.float32)

    for i, run in enumerate(runs):
        runid = run['runid']
        runkey = run['runkey']
        for ib, basin in enumerate(BASINS):
            mean_count[(runkey, basin)] = np.ma.asarray(annual_count_models[runkey][basin].mean(), dtype = np.float32)
            std_count[(runkey, basin)] = np.ma.asarray(annual_count_models[runkey][basin].std(), dtype = np.float32)
            for icat in np.arange(6):        
                mean_count_cat[(runkey, basin, icat)] = np.ma.asarray(annual_count_models[runkey][(basin,icat)].mean(), dtype = np.float32)
                std_count_cat[(runkey, basin, icat)] = np.ma.asarray(annual_count_models[runkey][(basin,icat)].std(), dtype = np.float32)

# do ACE
    for i, run in enumerate(runs):
        runid = run['runid']
        runkey = run['runkey']
        for ib, basin in enumerate(BASINS):
            #print('calc mean ace ',runid, runkey, basin)
            #print('annual_count no ace ', annual_count_models[runkey][('ace',basin)])

            mean_count[(runkey, 'ace', basin)] = np.ma.asarray(annual_count_models[runkey][('ace', basin)].mean(), dtype = np.float32)
            std_count[(runkey, 'ace', basin)] = np.ma.asarray(annual_count_models[runkey][('ace', basin)].std(), dtype = np.float32)
            for icat in np.arange(6):        
                mean_count_cat[(runkey, 'ace', basin, icat)] = np.ma.asarray(annual_count_models[runkey][('ace', basin, icat)].mean(), dtype = np.float32)
                std_count_cat[(runkey, 'ace', basin, icat)] = np.ma.asarray(annual_count_models[runkey][('ace', basin, icat)].std(), dtype = np.float32)


def _create_interann_plot(runids, runkeys, resols, annual_count_obs, \
                          annual_count_model, plot_basedir, BASINS, obs_datasets, colour_runid, colour_obs, years, std_count, plot_hurdat, plot_all_years=True, show_plots=False, use_model_years=False, ace = False):
    '''
    plot interannual variability compared to HURDAT observations
    '''

    # Set base part of the filename string for saving plots
    #figbase = '_'.join([r['runid'] for r in runs])
    cntlid = runids[0]
    testid = '_'.join([r for r in runids[1:]])
    #figbase = '_'.join([testid, cntlid])
    #if len(runids) > 10:
    #    figbase = '_'.join([runids[-1], cntlid])
    #else:
    #    figbase = '_'.join([testid, cntlid])
    figbase = '_'.join(runids)
    figbase = '_'.join(runkeys)

    alpha = 0.2

# plot the interannual variability
    for i, basin in enumerate(BASINS):        
        fig = plt.figure(figsize=(12,6),dpi=100)#dpi=300
        subpl=fig.add_subplot(111)
        lineplot_setup(annual_count_obs['years'], [0,BASINS[basin][3]], BASINS[basin][1],'Nominal Year','No. of TCs \n\n')
        fig = plt.gcf()

        minyear = 3000
        maxyear = 1000
        #print('run keys ',runkeys)
        for imod, runid in enumerate(runids):
            runkey = runkeys[imod]
            elts = runkey.split('_')
            algo = '_'.join(elts[1:])
            resol = resols[imod]
            if len(annual_count_model[runkey][basin]) > 12:
                if ace:
                    smoothed_data = tc_assess_code.smooth(np.array(annual_count_model[runkey]['ace',basin]), window_len=11, window='hanning')
                    std = std_count[(runkey, 'ace', basin)]
                else:
                    smoothed_data = tc_assess_code.smooth(np.array(annual_count_model[runkey][basin]), window_len=11, window='hanning')
                    std = std_count[(runkey, basin)]
                #print('nyears ',imod, runkey, len(annual_count_model[runkey][('years',BASINS[basin][0])]), annual_count_model[runkey][('years',BASINS[basin][0])])
                plt.plot(annual_count_model[runkey][('years',BASINS[basin][0])], smoothed_data, '-', \
                         label=''+runkey+'-'+resol+',%3.1f'%(std), color=colour_runid[imod], linewidth=3)

                if plot_all_years:
                    if ace:
                        plt.plot(annual_count_model[runkey][('years',BASINS[basin][0])], annual_count_model[runkey]['ace',basin], '-', color=colour_runid[imod], linewidth=3, alpha=alpha)
                    else:
                        plt.plot(annual_count_model[runkey][('years',BASINS[basin][0])], annual_count_model[runkey][basin], '-', color=colour_runid[imod], linewidth=3, alpha=alpha)
            else:
                if plot_all_years:
                    if ace:
                        plt.plot(annual_count_model[runkey][('years',BASINS[basin][0])], annual_count_model[runkey]['ace',basin], '-', color=colour_runid[imod], linewidth=3, \
                     alpha=alpha, label=''+runkey+'-'+resol+',%3.1f'%(std_count[(runkey,'ace',basin)]))
                    else:
                        #print('years for plot ',annual_count_model[runkey][('years',BASINS[basin][0])], annual_count_model[runkey][basin])
                        plt.plot(annual_count_model[runkey][('years',BASINS[basin][0])], annual_count_model[runkey][basin], '-', color=colour_runid[imod], linewidth=3, \
                     alpha=alpha, label=''+runkey+'-'+resol+',%3.1f'%(std_count[(runkey,basin)]))

            if len(annual_count_model[runkey][('years',BASINS[basin][0])]) > 1:
                minyear = min(min(annual_count_model[runkey][('years',BASINS[basin][0])]), minyear)
                maxyear = max(max(annual_count_model[runkey][('years',BASINS[basin][0])]), maxyear)
            else:
                minyear = annual_count_model[runkey][('years',BASINS[basin][0])][0]
                maxyear = annual_count_model[runkey][('years',BASINS[basin][0])][0]+1

        if plot_hurdat:
            window_len=11
            for io, obs_ds in enumerate(obs_datasets):
                if ace:
                    plt.plot(annual_count_obs['years'], annual_count_obs[obs_ds][('tc','ace',basin)], '-',color=colour_obs[io], linewidth=5, alpha=alpha)
                else:
                    plt.plot(annual_count_obs['years'], annual_count_obs[obs_ds][('tc',basin)], '-',color=colour_obs[io], linewidth=5, alpha=alpha)

                if len(annual_count_obs[obs_ds][('tc', basin)]) > window_len+1:
                    if ace:
                        smoothed_data = tc_assess_code.smooth(annual_count_obs[obs_ds][('tc', 'ace', basin)], window_len=window_len, window='hanning')
                    else:
                        smoothed_data = tc_assess_code.smooth(annual_count_obs[obs_ds][('tc', basin)], window_len=window_len, window='hanning')

                    plt.plot(annual_count_obs['years'], smoothed_data, '-',color=colour_obs[io], linewidth=5, label=obs_ds+'+') #+\
                 #',%3.1f'%(std_count[(obs_ds ,'all',basin)]))

                if ace:
                    plt.plot(annual_count_obs['years'], annual_count_obs[obs_ds][('hur','ace',basin)], '-.',color=colour_obs[io], linewidth=5, alpha=alpha)
                else:
                    plt.plot(annual_count_obs['years'], annual_count_obs[obs_ds][('hur',basin)], '-.',color=colour_obs[io], linewidth=5, alpha=alpha)

                if len(annual_count_obs[obs_ds][('hur', basin)]) > window_len+1:
                    if ace:
                        smoothed_data = tc_assess_code.smooth(annual_count_obs[obs_ds][('hur', 'ace', basin)], window_len=window_len, window='hanning')
                    else:
                        smoothed_data = tc_assess_code.smooth(annual_count_obs[obs_ds][('hur', basin)], window_len=window_len, window='hanning')

                    plt.plot(annual_count_obs['years'], smoothed_data, '-.',color=colour_obs[io], linewidth=5, label=obs_ds) #+\
                 #',%3.1f'%(std_count[(obs_ds, 'hur', basin)]))
        
        plt.legend(bbox_to_anchor=(1.0, 1.), loc=2, borderaxespad=-0.0, fontsize = 'small')
        if use_model_years:
            minyear = max(minyear, min(years))
            maxyear = min(maxyear, max(years))
        #print 'xlim ',minyear, maxyear
        subpl.set_xlim(minyear, maxyear)

        if ace:
            subpl.set_ylim(BASINS[basin][4], BASINS[basin][5])
            subpl.set_ylabel('TC ACE \n\n')
            plot_filename = os.path.join(plot_basedir, figbase+'_interann_var_'+algo+'_'+basin+'_ace.png')
        else:
            subpl.set_ylim(BASINS[basin][2], BASINS[basin][3])
            subpl.set_ylabel('TC frequency \n\n')
            plot_filename = os.path.join(plot_basedir, figbase+'_interann_var_'+algo+'_'+basin+'.png')
        fig.subplots_adjust(bottom=0.13, left=0.1, right=0.77, top=0.9, wspace=0.2, hspace=0.24)

        plt.savefig(plot_filename,dpi=100)
        #print 'save fig ',plot_filename
        if show_plots:
            plt.show()    

def calculate_correlation_years(years_model, data_model, years_obs, data_obs):
    #print years_model[0], years_model[-1]
    if years_model[0] in years_obs or years_model[-1] in years_obs and len(years_model) > 5:
        if years_model[0] >= years_obs[0]:
            # index in hurdat array for first year of model data
            index_start = np.where(years_obs == years_model[0])[0]
            index_end = np.where(years_obs == years_model[-1])[0]
            #print ('index_end ',index_start, index_end)
            #print('years obs for corr ',years_obs)
            #print('years obs for corr ',years_model)
            #print('data_model ',data_model)
            # if index_end has a value, then hurdat data is a superset of model data
            if len(index_end) == 1:
                hurdat_subset = data_obs[index_start[0]:index_end[0]+1]
                model_subset = data_model[0:len(hurdat_subset)]
                years_obs_subset = years_obs[index_start[0]:index_end[0]+1]
                years_model_subset = years_model[0:len(hurdat_subset)]
                hurdat_masked = np.ma.masked_where(np.ma.getmask(model_subset), hurdat_subset)
                if len(hurdat_masked) == len(model_subset):
                    xcorr = pearsonr(hurdat_masked, model_subset)
                else:
                    xcorr = np.array([0.0, 0.0])
                #print years_model_subset
                #print years_obs_subset
                #print model_subset
                #print hurdat_masked
                #print xcorr[0]
            else:
                # the last year of hurdat occurs before the last year of model
                index_model_end = np.where(years_obs[-1] == years_model)[0]
                #print 'index_model_end ',index_model_end[0]
                model_subset = data_model[0:index_model_end[0]+1]
                hurdat_subset = data_obs[index_start[0]:]
                hurdat_masked = np.ma.masked_where(np.ma.getmask(model_subset), hurdat_subset)
                if len(hurdat_masked) == len(model_subset):
                    xcorr = pearsonr(hurdat_masked, model_subset) 
                else:
                    xcorr = np.array([0.0, 0.0])
        else:
            # years_model[0] < years_hurdat[0]
            index_start = np.where(years_model == years_obs[0])[0]
            index_end = np.where(years_model == years_obs[-1])[0]
            # if index_end has a value, then hurdat data is a superset of model data
            if len(index_end) == 1:
                model_subset = data_model[index_start[0]:index_end[0]+1]
                hurdat_subset = data_obs
                hurdat_masked = np.ma.masked_where(np.ma.getmask(model_subset), hurdat_subset)
                if len(hurdat_masked) == len(model_subset):
                    xcorr = pearsonr(hurdat_masked, model_subset)
                else:
                    xcorr = np.array([0.0, 0.0])
            else:
                # the last year of hurdat occurs before the last year of model
                #print('index_end ',index_end)
                #print('years_model ',years_model)
                #print('years_obs ',years_obs)
                index_hurdat_end = np.where(years_model[-1] == years_obs)[0]
                #print('index_hurdat_end ',index_hurdat_end) 
                if len(index_hurdat_end) > 0:
                    hurdat_subset = data_obs[0:index_hurdat_end[0]+1]
                    model_subset = data_model[index_start[0]:]
                    hurdat_masked = np.ma.masked_where(np.ma.getmask(model_subset), hurdat_subset)
                    xcorr = pearsonr(hurdat_masked, model_subset) 
                else:
                    xcorr = np.array([0.0, 1.0])
    else:
        #print('model years not contained in hurdat years ',years_model[0])
        xcorr = np.array([0.0, 1.0])
    
    return xcorr

def _basin_correlations(runids, runkeys, BASINS, plot_basedir, annual_count_models, annual_count_obs, annual_corr_models):
    '''
    Calculate linear correlation coefficients for each model and each basin
    '''
    # Set base part of the filename string for saving plots
    #figbase = '_'.join([r['runid'] for r in runs])
    cntlid = runids[0]
    testid = '_'.join([r for r in runids[1:]])
    #figbase = '_'.join([testid, cntlid])
    figbase = '_'.join(runids)
    figbase = '_'.join(runkeys)

    for imod, runid in enumerate(runids):
        runkey = runkeys[imod]
        annual_corr_models[runkey] = {}

# first find out the overlapping years between models and observations - NH and SH may be slightly different
    for i, basin in enumerate(BASINS):        
        for imod, run in enumerate(runids):
            runkey = runkeys[imod]
            #resol = run['resol']
            years_model = annual_count_models[runkey][('years',BASINS[basin][0])]
            for metric in ['count','ace']:
                if metric == 'ace':
                    data_model = annual_count_models[runkey][('ace', basin)]
                else:
                    data_model = annual_count_models[runkey][basin]
            # indx = np.in1d(years_model, years_hurdat) # find [True, False] array where years equal
            
                obs_ds = 'hurdat2'
                for stype in ['hur','tc']:
                    years_obs = annual_count_obs['years']
                    if metric == 'ace':
                        data_obs = annual_count_obs[obs_ds][(stype,'ace',basin)]
                    else:
                        data_obs = annual_count_obs[obs_ds][(stype,basin)]

                    if len(years_model) > 1:
                        corr = calculate_correlation_years(years_model, data_model, years_obs, data_obs)
                    else:
                        print('correlation set to zero as years_model = 0', basin, run, runkey, metric, years_model)
                        corr = np.array([0.0, 1.0])
                    #print 'corr ', runkey, stype, basin, corr[0]
                
                    if metric == 'ace':
                        annual_corr_models[runkey].update({(basin,stype,'ace'): corr})
                    else:
                        annual_corr_models[runkey].update({(basin,stype): corr})
                
    for ib, basin in enumerate(BASINS):
        for imod, runid in enumerate(runids):
            runkey = runkeys[imod]
            filename = os.path.join(plot_basedir, runkey+'_interannual_'+basin+'_correlation.txt')    
            f = open(filename,'w')
            value = basin + ' '+runkey + ' '+str(annual_corr_models[runkey][(basin, 'tc')][0])+' ' + \
                    str(annual_corr_models[runkey][(basin, 'tc')][1])+' tc \n'
            f.write(value)
            value = basin + ' '+runkey + ' '+str(annual_corr_models[runkey][(basin, 'hur')][0])+' ' + \
                    str(annual_corr_models[runkey][(basin, 'hur')][1])+' hur \n'
            f.write(value)
            value = basin + ' '+runkey + ' '+str(annual_corr_models[runkey][(basin, 'tc', 'ace')][0])+' ' + \
                    str(annual_corr_models[runkey][(basin, 'tc', 'ace')][1])+' tc ace \n'
            f.write(value)
            value = basin + ' '+runkey + ' '+str(annual_corr_models[runkey][(basin, 'hur', 'ace')][0])+' ' + \
                    str(annual_corr_models[runkey][(basin, 'hur', 'ace')][1])+' hur ace \n'
            f.write(value)
            f.close()
    

def _create_correlation_plot(runids, runkeys, annual_corr_models, BASINS, \
                             BASIN_NAME, BASINS_LIST, plot_basedir, colour_runid, resols, show_plots=False, ace=False, MAX_RUNS = 4):
    '''
    create a histogram plot of the basin TC frequencies compared to HURDAT observations
    '''
    cntlid = runids[0]
    testid = '_'.join([r for r in runids[1:]])
    # Set base part of the filename string for saving plots
    figbase = '_'.join([testid, cntlid])
    if len(runids) > 10:
        figbase = '_'.join([runids[-1], cntlid])
    else:
        figbase = '_'.join([testid, cntlid])
    figbase = '_'.join(runids)
    figbase = '_'.join(runkeys)

    #print('basin_name ',BASIN_NAME)
    #print('basins_list ',BASINS_LIST)
    #print('basins ',BASINS)

    ylim = [-0.5,1]

    nbars = 3 + len(runids)  # hurdat*2+ctl+number of other runs
#    width = 0.1       # the width of the bars, wider when less runs
    width = 1.5*max(0.05, 1./float(nbars+0.5))
#    width = max(0.05, 1./float(nbars+3.5))
    
    if len(runids) <= MAX_RUNS:
        
        fig = plt.figure(figsize=(11,7),dpi=100)
        ax = fig.add_subplot(1,1,1)
        ind = np.arange(len(BASINS_LIST)) + 1
        for i, runid in enumerate(runids):
            resol = resols[i]
            runkey = runkeys[i]
            elts = runkey.split('_')
            algo = '_'.join(elts[1:])
            if ace:
                basin_corr_tc = [annual_corr_models[runkey][(basin, 'tc','ace')][0] for ib1, basin in enumerate(BASINS_LIST)]
                basin_corr_hur = [annual_corr_models[runkey][(basin, 'hur','ace')][0] for ib1, basin in enumerate(BASINS_LIST)]
            else:
                basin_corr_tc = [annual_corr_models[runkey][(basin, 'tc')][0] for ib1, basin in enumerate(BASINS_LIST)]
                basin_corr_hur = [annual_corr_models[runkey][(basin, 'hur')][0] for ib1, basin in enumerate(BASINS_LIST)]
            #print 'basin_corr_tc ',basin_corr_tc
            
            ax.barh(ind+width*i, basin_corr_tc, width/3., \
                    color=colour_runid[i], label=''+runkey+'-'+resol, alpha = 0.5)
            ax.barh(ind+width*(i+0.5), basin_corr_hur, width/3., \
                    color=colour_runid[i], label='HUR')
        
            plt.xlim(-0.5, 1)
            ax.set_yticks(ind+width*nbars/2)
            ax.set_yticklabels( (BASINS_LIST) )
            ax.set_ylim(ind.min()-1, ind.max()+1)
            plt.gca().invert_yaxis()
            if ace:
                ax.set_xlabel('ACE correlation \n\n')
                plt.title('TC ACE correlation')
            else:
                ax.set_xlabel('Frequency correlation \n\n')
                plt.title('TC frequency correlation')

            plt.legend(bbox_to_anchor=(1.0, 1.), loc=2, borderaxespad=-0.0, fontsize='medium')
        fig.subplots_adjust(bottom=0.1, left=0.15, right=0.74, top=0.94, wspace=0.2, hspace=0.24)
            
        if ace:
            plot_filename = os.path.join(plot_basedir, figbase+'_basin_correlation_ace_'+algo+'.png')
        else:
            plot_filename = os.path.join(plot_basedir, figbase+'_basin_correlation_'+algo+'.png')
        plt.savefig(plot_filename)
        if show_plots:
            plt.show()    
        
    else:
        width = 0.4
        for ib, basin in enumerate(BASINS_LIST):
            fig = plt.figure(figsize=(10,9),dpi=100)
            ax = fig.add_subplot(1,1,1)
            ind = np.arange(len([basin]))
            ind_runid = np.arange(len(runids))
            run_labels = []; run_colour = []; run_number = []; run_resol =[]
            b_corr_tc = []; b_corr_hur = []
            for i, runid in enumerate(runids):
                runkey = runkeys[i]
                resol = resols[i]
                runids = [r for r in runids]
                if ace:
#                    basin_corr_tc = [annual_corr_models[runkey][(bas, 'tc','ace')][0] for ib1, bas in enumerate([basin])]
#                    basin_corr_hur = [annual_corr_models[runkey][(bas, 'hur','ace')][0] for ib1, bas in enumerate([basin])]
                    basin_corr_tc = annual_corr_models[runkey][(basin, 'tc','ace')][0]
                    basin_corr_hur = annual_corr_models[runkey][(basin, 'hur','ace')][0]
                else:
#                    basin_corr_tc = [annual_corr_models[runkey][(bas, 'tc')][0] for ib1, bas in enumerate([basin])]
#                    basin_corr_hur = [annual_corr_models[runkey][(bas, 'hur')][0] for ib1, bas in enumerate([basin])]
                    basin_corr_tc = annual_corr_models[runkey][(basin, 'tc')][0]
                    basin_corr_hur = annual_corr_models[runkey][(basin, 'hur')][0]
                b_corr_tc.append(basin_corr_tc)
                b_corr_hur.append(basin_corr_hur)
                #if run['ensemble_member'] != '':
                #    run_labels.append(run['ensemble_member'])
                #else:
                #    run_labels.append(run['runkey'])
                run_labels.append(runkey)

                run_colour.append(colour_runid[i])
                run_number.append(i)
                #print(runid)
                #print run
                run_resol.append(resol)
            #print len(run_colour), ind_runid
            #print 'b_corr_tc ',b_corr_tc

            ax.barh(ind_runid, np.asarray(b_corr_tc), width, color = run_colour, alpha=0.5)
            ax.barh(ind_runid+width, np.asarray(b_corr_hur), width, \
                       color=run_colour)
            
            #ax.set_yticks(ind+width*nbars/2)
            #ax.set_yticks(ind_runid)
            plt.yticks(ind_runid+width, run_number)
            #ax.set_yticklabels( (run_labels) )
            plt.xlim(-0.2, 1)
            ax.set_ylim(ind_runid.min(), ind_runid.max()+1)
            if ace:
                ax.set_xlabel('ACE correlation \n\n')
            else:
                ax.set_xlabel('Frequency correlation \n\n')
            ax2 = ax.twinx()
            plt.yticks(ind_runid+width)
            ax2.set_ylim(ind_runid.min(), ind_runid.max()+1)
            ax2.set_yticklabels(run_labels)
            #ax2.set_yticks(ind_runid)
            #plt.yticks(ind_runid+width, run_labels, fontsize='small')
            #plt.gca().invert_yaxis()
            plt.plot([-0.2,1], [ind_runid.max()-2, ind_runid.max()-2], 'black')
            plt.legend(bbox_to_anchor=(1.0, 1.), loc=2, borderaxespad=-0.0, fontsize='small')
            if ace:
                plt.title('TC ACE correlation - '+BASIN_NAME[ib]+'('+basin+')')
            else:
                plt.title('TC frequency correlation - '+BASIN_NAME[ib]+'('+basin+')')
            fig.subplots_adjust(bottom=0.1, left=0.1, right=0.73, top=0.94, wspace=0.2, hspace=0.24)
            
            if ace:
                plot_filename = os.path.join(plot_basedir, figbase+'_'+basin+'_correlation_ace_'+algo+'.png')
            else:
                plot_filename = os.path.join(plot_basedir, figbase+'_'+basin+'_correlation_'+algo+'.png')
            plt.savefig(plot_filename)
            if show_plots:
                plt.show()    
        
        
    #plt.close()

def _create_basin_plot(runids, runkeys, mean_count, annual_count_obs, BASINS, \
                             BASIN_NAME, BASINS_LIST, plot_basedir, colour_runid, resols, show_plots=False, ace=False, MAX_RUNS = 4):
    '''
    create a histogram plot of the basin TC frequencies compared to HURDAT observations
    '''
    cntlid = runids[0]
    testid = '_'.join([r for r in runids[1:]])
    # Set base part of the filename string for saving plots
    figbase = '_'.join([testid, cntlid])
    if len(runids) > 10:
        figbase = '_'.join([runids[-1], cntlid])
    else:
        figbase = '_'.join([testid, cntlid])
    figbase = '_'.join(runids)
    figbase = '_'.join(runkeys)
 
    ylim = [-0.5,1]

    nbars = 3 + len(runids)  # hurdat*2+ctl+number of other runs
#    width = 0.1       # the width of the bars, wider when less runs
    width = 1.*max(0.05, 1./float(nbars+0.5))
#    width = max(0.05, 1./float(nbars+3.5))
    
    if len(runids) <= MAX_RUNS:
        
        fig = plt.figure(figsize=(11,7),dpi=100)
        ax = fig.add_subplot(1,1,1)
        ind = np.arange(len(BASINS_LIST))+1.0

        basin_count_obs = []
        for ib1, basin in enumerate(BASINS_LIST):
            #print('obs val ',annual_count_obs['hurdat2'][('tc', basin)])
            count = np.ma.asarray(annual_count_obs['hurdat2'][('tc', basin)].mean(), dtype = np.float32)
            basin_count_obs.append(count)
        
        ax.barh(ind, basin_count_obs, width/2., \
                color='black', label='OBS', alpha = 0.8)

        for i, runid in enumerate(runids):
            resol = resols[i]
            runkey = runkeys[i]
            elts = runkey.split('_')
            algo = '_'.join(elts[1:])
            #print('runkey ',runkeys, algo)
            basin_count = []
            for ib1, basin in enumerate(BASINS_LIST):
                if ace:
                    count = mean_count[(runkey, 'ace', basin)]
                    basin_count.append(count)
                else:
                    basin_count.append(mean_count[(runkey,basin)])
            #print('basin_count ',runid, basin, basin_count)
            #print 'basin_corr_tc ',basin_corr_tc
            
            ax.barh(ind+width*(i+1), basin_count, width/2., \
                    color=colour_runid[i], label=''+runkey+'-'+resol, alpha = 0.8)
            plt.xlim(0, 60)
            ax.set_yticks(ind+width*nbars/3)
            ax.set_yticklabels( (BASINS_LIST) )
            ax.set_ylim(ind.min()-1, ind.max()+1)
            plt.gca().invert_yaxis()
            if ace:
                ax.set_xlabel('ACE count \n\n')
                plt.title('TC ACE count')
            else:
                ax.set_xlabel('Frequency count \n\n')
                plt.title('TC frequency count')

            plt.legend(bbox_to_anchor=(1.0, 1.), loc=2, borderaxespad=-0.0, fontsize='medium')
        fig.subplots_adjust(bottom=0.1, left=0.15, right=0.74, top=0.94, wspace=0.2, hspace=0.24)
            
        if ace:
            plot_filename = os.path.join(plot_basedir, figbase+'_basin_frequency_ace_'+algo+'.png')
        else:
            plot_filename = os.path.join(plot_basedir, figbase+'_basin_frequency_'+algo+'.png')
        plt.savefig(plot_filename)
        if show_plots:
            plt.show()    
        
    else:
        width = 0.4
        for ib, basin in enumerate(BASINS_LIST):
            fig = plt.figure(figsize=(10,9),dpi=100)
            ax = fig.add_subplot(1,1,1)
            ind = np.arange(len([basin]))
            ind_runid = np.arange(len(runids))
            run_labels = []; run_colour = []; run_number = []; run_resol =[]
            b_corr_tc = []; b_corr_hur = []
            for i, runid in enumerate(runids):
                runkey = runkeys[i]
                resol = resols[i]
                runids = [r for r in runids]
                if ace:
#                    basin_corr_tc = [annual_corr_models[runkey][(bas, 'tc','ace')][0] for ib1, bas in enumerate([basin])]
#                    basin_corr_hur = [annual_corr_models[runkey][(bas, 'hur','ace')][0] for ib1, bas in enumerate([basin])]
                    basin_corr_tc = annual_corr_models[runkey][(basin, 'tc','ace')][0]
                    basin_corr_hur = annual_corr_models[runkey][(basin, 'hur','ace')][0]
                else:
#                    basin_corr_tc = [annual_corr_models[runkey][(bas, 'tc')][0] for ib1, bas in enumerate([basin])]
#                    basin_corr_hur = [annual_corr_models[runkey][(bas, 'hur')][0] for ib1, bas in enumerate([basin])]
                    basin_corr_tc = annual_corr_models[runkey][(basin, 'tc')][0]
                    basin_corr_hur = annual_corr_models[runkey][(basin, 'hur')][0]
                b_corr_tc.append(basin_corr_tc)
                b_corr_hur.append(basin_corr_hur)
                #if run['ensemble_member'] != '':
                #    run_labels.append(run['ensemble_member'])
                #else:
                #    run_labels.append(run['runkey'])
                run_labels.append(runkey)

                run_colour.append(colour_runid[i])
                run_number.append(i)
                #print(runid)
                #print run
                run_resol.append(resol)
            #print len(run_colour), ind_runid
            #print 'b_corr_tc ',b_corr_tc

            ax.barh(ind_runid, np.asarray(b_corr_tc), width, color = run_colour, alpha=0.5)
            ax.barh(ind_runid+width, np.asarray(b_corr_hur), width, \
                       color=run_colour)
            
            #ax.set_yticks(ind+width*nbars/2)
            #ax.set_yticks(ind_runid)
            plt.yticks(ind_runid+width, run_number)
            #ax.set_yticklabels( (run_labels) )
            plt.xlim(-0.2, 1)
            ax.set_ylim(ind_runid.min(), ind_runid.max()+1)
            if ace:
                ax.set_xlabel('ACE correlation \n\n')
            else:
                ax.set_xlabel('Frequency correlation \n\n')
            ax2 = ax.twinx()
            plt.yticks(ind_runid+width)
            ax2.set_ylim(ind_runid.min(), ind_runid.max()+1)
            ax2.set_yticklabels(run_labels)
            #ax2.set_yticks(ind_runid)
            #plt.yticks(ind_runid+width, run_labels, fontsize='small')
            #plt.gca().invert_yaxis()
            plt.plot([-0.2,1], [ind_runid.max()-2, ind_runid.max()-2], 'black')
            plt.legend(bbox_to_anchor=(1.0, 1.), loc=2, borderaxespad=-0.0, fontsize='small')
            if ace:
                plt.title('TC ACE correlation - '+BASIN_NAME[ib]+'('+basin+')')
            else:
                plt.title('TC frequency correlation - '+BASIN_NAME[ib]+'('+basin+')')
            fig.subplots_adjust(bottom=0.1, left=0.1, right=0.73, top=0.94, wspace=0.2, hspace=0.24)
            
            if ace:
                plot_filename = os.path.join(plot_basedir, figbase+'_'+basin+'_frequency_ace.png')
            else:
                plot_filename = os.path.join(plot_basedir, figbase+'_'+basin+'_frequency.png')
            plt.savefig(plot_filename)
            if show_plots:
                plt.show()    
        

def tc_driver(runs, years_hurdat, years_model, plot_basedir, HURDAT_DATA_BASEDIR, do_metric=False, show_plots = False, use_model_years=False, plot_hurdat = True, plot_all_years = True, paired = False, do_new = False, do_new_obs = False):
    '''
    Routine to calculate and plot basin metrics of tropical cyclones.

    This is a type 2 multi-function routine that returns no metrics

    Arguments:
        runs - list of run dictionaries.  Each dictionary contains
               metadata for a single model run.  The first dictionary
               in this list is the control experiment.
               (see assessment_area.Area.py_program for description of
               the contents of this dictionary)

    Returns:
        if do_metric is True, returns a dictionary of metrics for the Autoassess metric plot
        if do_metric is false:
            doesn't return any objects - it only writes image files to the
            current working dir
    '''
    time_ref = time.perf_counter()

    return_code = 0
    colour_runid = {}
    for i, run in enumerate(runs):
        track_method = ''
        try:
            track_method = run['track_method']
        except:
            pass
        #runkey = run['runid']+run['ensemble_member']+track_method
        # run['area'] = runid+resol
        runkey = run['area']+run['ensemble_member']+track_method
        run.update({'runkey': runkey})
        if paired:
            colour_runid[runkey] = COLOURNAME_paired[i]
        else:
            colour_runid[runkey] = COLOURNAME_unpaired[i]
        #print i, run, runkey, colour_runid[runkey]
    years = years_model

    # initialise dictionaries to hold the mean count over the period, and the standard deviation
    mean_count = {}
    std_count = {}
    mean_count_cat = {}
    std_count_cat = {}

    annual_count_models = {}
    annual_corr_models = {}
    density_cube = {}
    genesis_cube = {}

    annual_count_obs = {}
    density_cube_obs = {}
    genesis_cube_obs = {}
    
    annual_count_obs['years'] = years_hurdat

    for obs_ds in OBS_DATASETS:
        annual_count_obs[obs_ds] = {}
        genesis_cube_obs[obs_ds] = {}
        density_cube_obs[obs_ds] = {}
        filename_obs_all_and_hur = obs_ds+'_alltcs_hur_%s_%s.pkl' % (years_hurdat[0], years_hurdat[-1])
        #print 'obs_path ',os.path.join(HURDAT_DATA_BASEDIR, filename_obs_all_and_hur)

# calculate the number of storms in each basin for HURDAT
        if not (os.path.exists(os.path.join(HURDAT_DATA_BASEDIR, filename_obs_all_and_hur))) or do_new_obs:
            if obs_ds == 'hurdat2':
                obs_storms = ts_obs.load_data.load_data()
            else:
                obs_storms = ts_obs.load_data_ibtracs.load_data()

            calculate_hurdat_frequency(obs_storms, annual_count_obs, obs_ds)
            
# calculate the track densities for each hemisphere
            track_density = None; genesis_density = None
            for type, cube_dict in zip(['genesis','track'],[genesis_cube_obs[obs_ds], density_cube_obs[obs_ds]]):
                if type == 'genesis': do_genesis = True
                else: do_genesis = False
                calculate_hurdat_density(obs_storms, annual_count_obs, cube_dict, obs_ds, genesis=do_genesis)

            fh = open(os.path.join(HURDAT_DATA_BASEDIR,filename_obs_all_and_hur), 'wb') # write binary file mode
            pickle.dump(annual_count_obs[obs_ds], fh)
            pickle.dump(density_cube_obs[obs_ds], fh)
            pickle.dump(genesis_cube_obs[obs_ds], fh)
            fh.close()
    
        else:
    #reading
            #print ('obs file ',os.path.join(HURDAT_DATA_BASEDIR,filename_obs_all_and_hur))
            fh = open(os.path.join(HURDAT_DATA_BASEDIR,filename_obs_all_and_hur), 'rb') #read binary mode
            annual_count_obs[obs_ds] = pickle.load(fh)
            density_cube_obs[obs_ds] = pickle.load(fh)
            genesis_cube_obs[obs_ds] = pickle.load(fh)
            fh.close()
    #print 'obs density cube ',density_cube_obs[obs_ds]
    
# loop over 2 types of hurdat, all and hurr
    for obs_ds in OBS_DATASETS:
        for i, key in enumerate(OBS_TYPES):
            density_cube[(obs_ds, key)] = density_cube_obs[obs_ds][key] 
        for i, key in enumerate(OBS_TYPES):
            genesis_cube[(obs_ds, key)] = genesis_cube_obs[obs_ds][key] 
    #print 'density_cube_hurdat',density_cube

    time_obs = time.perf_counter()
    print('Var Reading obs took ',(time_obs - time_ref),' seconds')

# calculate the number of storms in each basin for each model
    for i, run in enumerate(runs):
        runkey = run['runkey']
        runid, filename_pickle, fnames_model = discover_model_track_files(run, annual_count_models, years_model, use_model_years)
        print(runid, runkey)

        nensemble = run['ens_number']
        #print ('pkl file ', os.path.join(run['data_root'], run['area'], filename_pickle))
        if not os.path.exists(os.path.join(run['data_root'], run['area'], filename_pickle)) or do_new:
# if want clean run, then ignore this
            storms_model_all = load_example_storms(fnames_model, run['wind10m'])
            if use_model_years:
                storms_model_all = storm_in_years(storms_model_all, years_model)
            storms_model_cat0 = filter_storms(storms_model_all, sfilter = [0])
            storms_model = storms_model_all
            calculate_model_frequency(fnames_model, runid, runkey, annual_count_models, nensemble, storms_model)
            calculate_model_ace(fnames_model, runid, runkey, annual_count_models, nensemble, storms_model)

            fh = open(os.path.join(run['data_root'], run['area'], filename_pickle), 'wb') # write binary file mode
            pickle.dump(annual_count_models[runkey], fh)
            pickle.dump(storms_model, fh)
            pickle.dump(storms_model_cat0, fh)
            fh.close()
        else:
            fh = open(os.path.join(run['data_root'], run['area'], filename_pickle), 'rb') #read binary mode
            annual_count_models[runkey] = pickle.load(fh)
            storms_model = pickle.load(fh)
            storms_model_cat0 = pickle.load(fh)
            fh.close()

    time_model = time.perf_counter()
    print('Var Model read took ',(time_model - time_obs),' seconds')

    calc_tc_frequency_mean_std(runs, BASINS, annual_count_obs, OBS_DATASETS, annual_count_models, OBS_TYPES, mean_count, std_count, \
                               mean_count_cat, std_count_cat)

    time_model_mean = time.perf_counter()
    print('Var Model mean took ',(time_model_mean - time_model),' seconds')

   # interannual metrics calculation, correlation etc (only if AMIP, how do we know if it AMIP?)
    metrics = _basin_frequency_metrics(runs, mean_count, std_count, mean_count_cat, \
                                std_count_cat, BASINS, BASIN_NAME[:len(BASINS)], plot_basedir, annual_count_models, annual_count_obs)

    time_metrics = time.perf_counter()
    print('Var Model mean took ',(time_metrics - time_model_mean),' seconds')

    # calculate correlations
    correlations = _basin_correlations(runs, BASINS, BASIN_NAME, plot_basedir, annual_count_models, annual_count_obs, \
                                       annual_corr_models)

    time_corr = time.perf_counter()
    print('Var Model corr took ',(time_corr - time_metrics),' seconds')

    if not do_metric:
        _create_correlation_plot(runs, annual_corr_models, BASINS, \
                       BASIN_NAME, plot_basedir, colour_runid, show_plots)

        _create_correlation_plot(runs, annual_corr_models, BASINS, \
                       BASIN_NAME, plot_basedir, colour_runid, show_plots, ace=True)
# basin histogram plot
    if not do_metric:
        _create_basin_plot(runs, mean_count, std_count, mean_count_cat, std_count_cat, BASINS, BASIN_NAME[:len(BASINS)], \
                           plot_basedir, colour_runid, show_plots, paired = paired)

        _create_basin_plot(runs, mean_count, std_count, mean_count_cat, std_count_cat, BASINS, BASIN_NAME[:len(BASINS)], \
                           plot_basedir, colour_runid, show_plots, paired = paired, ace = True)

        #_create_basin_plot_horiz(runs, mean_count, std_count, mean_count_cat, std_count_cat, BASINS, BASIN_NAME[:len(BASINS)], \
        #                   plot_basedir, colour_runid, show_plots, plot_cat0=True)

# interannual variability plot
    if not do_metric:
        _create_interann_plot(runs, annual_count_obs, \
                    annual_count_models, plot_basedir, colour_runid, years, std_count, show_plots, plot_hurdat, plot_all_years, use_model_years=use_model_years)
        _create_interann_plot(runs, annual_count_obs, \
                    annual_count_models, plot_basedir, colour_runid, years, std_count, show_plots, plot_hurdat, plot_all_years, use_model_years=use_model_years, ace = True)
    
# track density plot
    if not do_metric:
        _create_density_ncfiles(runs, density_cube, plot_basedir, show_plots)

        _create_genesis_ncfiles(runs, genesis_cube, plot_basedir, show_plots)

        for dentype, cube in zip(['den','gen'], [density_cube, genesis_cube]):
            for plot_type in ['', 'ts','diff']:
                _create_density_plot(runs, cube, plot_basedir, show_plots, plot_type = plot_type, title = dentype)

    time_plots = time.perf_counter()
    print('Var Model corr took ',(time_plots - time_corr),' seconds')

    return metrics

def main_metric(run):
    #print type(run)
    print(run)
    # easier to have array to enable enumerate across dictionary(ies)
    runs = [run]
#    plot_basedir = '/data/local/hadom/maverick/plots/'
    plot_basedir = os.getcwd()
    # want hurdat years to encompass model years but can be longer, e.g. standard 20-30yr period
    years_hurdat = np.arange(1979,2011)
    years_model = np.arange(1979,2079)
    
    #print 'runs ',runs
    #print run['runid']

    metrics = tc_driver(runs, years_hurdat, years_model, plot_basedir, HURDAT_DATA_BASEDIR, do_metric=True)
    #print 'metrics ',metrics
    return metrics
    
def main_comparison(runs):
#    plot_basedir = '/data/local/hadom/maverick/plots/'
    plot_basedir = os.getcwd()
    # want hurdat years to encompass model years but can be longer, e.g. standard 20-30yr period
    years_hurdat = np.arange(1979,2011)
    years_model = np.arange(1979,2011)
    
    #print 'runs ',runs
    #print run['runid']
    show_plots = False

    metrics = tc_driver(runs, years_hurdat, years_model, plot_basedir, HURDAT_DATA_BASEDIR, show_plots=show_plots)
    print('metrics returned ',metrics)

if __name__ == '__main__':
    plot_basedir = '/data/local/hadom/maverick/plots/'
    # want hurdat years to encompass model years but can be longer, e.g. standard 20-30yr period
    years_hurdat = np.arange(1979,2011)
    years_model = np.arange(1979,2011)

#    return_code = tc_driver(runs, years_hurdat, years_model, plot_basedir, HURDAT_DATA_BASEDIR)
    return_code = main(runs)
    
    print(return_code)
