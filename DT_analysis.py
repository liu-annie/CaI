##
drive_path = 'c:/'
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from scipy.stats import anderson_ksamp
from scipy.stats import kruskal
from scipy.stats import variation
from scipy import signal as sps
import seaborn as sns
import glob
import re
##
#This piece spits out all the peaks in one dataframe
def getpeaks(date):
    '''Spits out all the peaks from imaging session
    session input as string
    '''
    # This piece spits out all the peaks from one session in one dataframe
    peakdf = pd.DataFrame([])
    os.chdir('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\HabituationFiles\\%s' % date)
    for filename in glob.glob('*dt.txt'):
        f = pd.read_csv(filename, nrows=175)
        df = f[[col for col in f.columns if 'G PMT' in col]]
        peak = []
        for col in df.columns:
            a = df[col]
            firsta = 1;
            firstb = 24;
            # Figures out if there is a min or max and sees if it passes threshold (3SD)
            if np.absolute(min(a[26:80])) > np.absolute(max(a[26:80])) and np.absolute(min(a[26:80])) >= 3 * np.std(
                    df[col][firsta:firstb]):
                b = min(a[26:80])
                peak.append(b)
            elif np.absolute(max(a[26:80])) > np.absolute(min(a[26:80])) and np.absolute(max(a[26:80])) >= 3 * np.std(
                    df[col][firsta:firstb]):
                b = max(a[26:80])
                peak.append(b)
            else:
                b = 0
                peak.append(b)
            peaks = pd.DataFrame(peak).T
        peaks.columns = df.columns
        peaks = pd.concat([pd.DataFrame({'Trial': [int(filename.split('dt')[0])]}), peaks], axis=1)
        peakdf = peakdf.append(peaks, ignore_index=True)
    peakdf.to_csv('%s_peaks.csv' % date, index=False)
    trials = pd.read_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\Analysis\\Odor_Panel\\Odor_Trials.csv')
    filerow = trials.loc[trials['File'] == date]
    odortrials = {}
    for t in filerow.Odor.unique():
        y = {t: [int(x) for x in filerow.loc[filerow['Odor'] == t][['T1', 'T2', 'T3', 'T4']].values.tolist()[0]]}
        odortrials.update(y)
    # Get average peak across all trials using peakdf dataframe
    meandf = pd.DataFrame([])
    for key in odortrials:
        odor = odortrials[key]
        mean = []
        for col in peakdf.loc[peakdf['Trial'].isin(odor)][
            [col for col in peakdf.loc[peakdf['Trial'].isin(odor)].columns if 'G PMT' in col]]:
            mean.append(peakdf.loc[peakdf['Trial'].isin(odor)][col].mean())
        mean = pd.DataFrame(mean).T
        mean.columns = peakdf.loc[peakdf['Trial'].isin(odor)][
            [col for col in peakdf.loc[peakdf['Trial'].isin(odor)].columns if 'G PMT' in col]].columns
        meandf = meandf.append(mean)
    meandf = meandf.reset_index(drop=True)
    meandf.columns = [str(col) + '_' + date for col in meandf.columns]
    meandf = pd.concat([pd.DataFrame({'Odor': odortrials.keys()}), meandf], axis=1)
    meandf.to_csv('%s_mean.csv' % date, index=False)
    # Get proportion of successful trials
    successdf = pd.DataFrame([])
    for key in odortrials:
        odor = odortrials[key]
        newdf = peakdf.loc[peakdf['Trial'].isin(odor)]
        s = []
        for col in peakdf.loc[peakdf['Trial'].isin(odor)][
            [col for col in peakdf.loc[peakdf['Trial'].isin(odor)].columns if 'G PMT' in col]]:
            s.append(np.divide((newdf.loc[:, col] != 0).sum(), float(len(newdf.loc[:, col]))))
        s = pd.DataFrame(s).T
        s.columns = peakdf.loc[peakdf['Trial'].isin(odor)][
            [col for col in peakdf.loc[peakdf['Trial'].isin(odor)].columns if 'G PMT' in col]].columns
        successdf = successdf.append(s)
    successdf = successdf.reset_index(drop=True)
    successdf.columns = [str(col) + '_' + date for col in successdf.columns]
    successdf = pd.concat([pd.DataFrame({'Odor': odortrials.keys()}), successdf], axis=1)
    successdf.to_csv('%s_success.csv' % date, index=False)
    return 'Done'
##
def getintegral(date):
    '''Compute integrals and integral means
    date: string, session
    '''
    temp = pd.DataFrame([])
    os.chdir('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\HabituationFiles\\%s' % date)
    # Pull the trials that correspond to specific date/odors
    trials = pd.read_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\Analysis\\Odor_Panel\\Odor_Trials.csv')
    filerow = trials.loc[trials['File'] == date]
    odortrials = {}
    for t in filerow.Odor.unique():
        y = {t: [int(x) for x in filerow.loc[filerow['Odor'] == t][['T1', 'T2', 'T3', 'T4']].values.tolist()[0]]}
        odortrials.update(y)
    # Get the frame rate for a specified date
    num = trials.File.unique().tolist().index('%s' % date)
    fr = trials.loc[trials['File'] == trials.File.unique().tolist()[num]]['FrameRate'].iloc[0]
    # Get the integral
    intdf = pd.DataFrame([])
    for filename in glob.glob('*dt.txt'):
        f = pd.read_csv(filename, nrows=125)
        df = f[[col for col in f.columns if 'G PMT' in col]]
        winstart = np.int(4 * fr)
        winend = np.int(12 * fr)
        integral = []
        for col in df.columns:
            a = df[col]
            firsta = 1;
            firstb = 24;
            # Figures out if there is a min or max and sees if it passes threshold (3SD)
            if np.absolute(min(a[26:80])) > np.absolute(max(a[26:80])) and np.absolute(min(a[26:80])) >= 3 * np.std(
                    df[col][firsta:firstb]):
                b = sum(df[col][winstart:winend] * (1 / fr))
                integral.append(b)
            elif np.absolute(max(a[26:80])) > np.absolute(min(a[26:80])) and np.absolute(max(a[26:80])) >= 3 * np.std(
                    df[col][firsta:firstb]):
                b = sum(df[col][winstart:winend] * (1 / fr))
                integral.append(b)
            else:
                b = 0
                integral.append(b)
        integral = pd.DataFrame(integral).T
        integral.columns = df.columns
        integral = pd.concat([pd.DataFrame({'Trial': [int(filename.split('dt')[0])]}), integral], axis=1)
        intdf = intdf.append(integral)
    intdf.to_csv('%s_integral.csv' % date, index=False)
    # Get average integral across all trials using integral dataframe
    meanint = pd.DataFrame([])
    for key in odortrials:
        odor = odortrials[key]
        mean = []
        for col in intdf.loc[intdf['Trial'].isin(odor)][
            [col for col in intdf.loc[intdf['Trial'].isin(odor)].columns if 'G PMT' in col]]:
            mean.append(intdf.loc[intdf['Trial'].isin(odor)][col].mean())
        mean = pd.DataFrame(mean).T
        mean.columns = intdf.loc[intdf['Trial'].isin(odor)][
            [col for col in intdf.loc[intdf['Trial'].isin(odor)].columns if 'G PMT' in col]].columns
        meanint = meanint.append(mean)
    meanint = meanint.reset_index(drop=True)
    meanint.columns = [str(col) + '_' + date for col in meanint.columns]
    meanint = pd.concat([pd.DataFrame({'Odor': odortrials.keys()}), meanint], axis=1)
    meanint.to_csv('%s_meanint.csv' % date, index=False)
    return 'Done'
##
def getbaseline(date):
    temp = pd.DataFrame([])
    os.chdir('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\HabituationFiles\\%s' % date)
    # Pull the trials that correspond to specific date/odors
    trials = pd.read_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\Analysis\\Odor_Panel\\Odor_Trials.csv')
    filerow = trials.loc[trials['File'] == date]
    odortrials = {}
    for t in filerow.Odor.unique():
        y = {t: [int(x) for x in filerow.loc[filerow['Odor'] == t][['T1', 'T2', 'T3', 'T4']].values.tolist()[0]]}
        odortrials.update(y)
    # Get the frame rate for a specified date
    num = trials.File.unique().tolist().index('%s' % date)
    fr = trials.loc[trials['File'] == trials.File.unique().tolist()[num]]['FrameRate'].iloc[0]
    # Get baseline
    baseline = pd.DataFrame([])
    for filename in glob.glob('*dt.txt'):
        f = pd.read_csv(filename, nrows=125)
        df = f[[col for col in f.columns if 'G PMT' in col]]
        winstart = np.int(4 * fr)
        winend = np.int(12 * fr)
        base = []
        for col in df.columns:
            a = df[col]
            firsta = 1;
            firstb = 24;
            b = (df[col][firsta:firstb]).mean()
            base.append(b)
        base = pd.DataFrame(base).T
        base.columns = df.columns
        base = pd.concat([pd.DataFrame({'Trial': [int(filename.split('dt')[0])]}), base], axis=1)
        baseline = baseline.append(base)
    baseline.to_csv('%s_baseline.csv' % date, index=False)
    # mean baseline
    meanbase = pd.DataFrame([])
    for key in odortrials:
        odor = odortrials[key]
        mean = []
        for col in baseline.loc[baseline['Trial'].isin(odor)][
            [col for col in baseline.loc[baseline['Trial'].isin(odor)].columns if 'G PMT' in col]]:
            mean.append(baseline.loc[baseline['Trial'].isin(odor)][col].mean())
        mean = pd.DataFrame(mean).T
        mean.columns = baseline.loc[baseline['Trial'].isin(odor)][
            [col for col in baseline.loc[baseline['Trial'].isin(odor)].columns if 'G PMT' in col]].columns
        meanbase = meanbase.append(mean)
    meanbase = meanbase.reset_index(drop=True)
    meanbase.columns = [str(col) + '_' + date for col in meanbase.columns]
    meanbase = pd.concat([pd.DataFrame({'Odor': odortrials.keys()}), meanbase], axis=1)
    meanbase.to_csv('%s_meanbase.csv'%date,index=False)
    return 'Done'
##
def concat(odorfile):
    '''trials is the file that has odor trial orders
    'C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\Analysis\\Odor_Panel\\Odor_Trials.csv'
    'C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\Analysis\\Odor_Panel\\Habituation_Trials.csv'
    '''
    trials = pd.read_csv(odorfile)
    filerow = trials.loc[trials['File'] == date]
    odortrials = {}
    for t in trials.Odor.unique():
        y = {t: [int(x) for x in filerow.loc[filerow['Odor'] == t][['T1', 'T2', 'T3', 'T4']].values.tolist()[0]]}
        odortrials.update(y)
    fullpeak = pd.DataFrame([])
    # Concat peak responses
    for date in trials.File.unique():
        # reorganize dataframes
        mean = pd.read_csv(
            'C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\{0}\\{1}_mean.csv'.format(unicode(date, 'utf-8'),
                                                                                                 unicode(date,
                                                                                                         'utf-8')))
        mdf = pd.concat([mean['Odor'], mean[[col for col in mean.columns if ')_' in col]]], axis=1)
        temp = mdf.T
        temp.reset_index(level=0, inplace=True)
        temp.columns = temp.iloc[0]
        temp = temp.reindex(temp.index.drop(0))
        temp.rename(columns={'Odor': 'Mouse'}, inplace=True)
        temp = temp.reset_index(drop=True)
        group = []
        for x in list(temp.index.values):
            temp.iloc[x]['Mouse'] = temp.iloc[x]['Mouse'].split(')_')[1]
            indexnum = trials.File.unique().tolist().index(temp['Mouse'][x])
            groupname = trials.loc[trials.File == trials.File.unique()[indexnum]].Group.iloc[0]
            group.append(groupname)
        group = pd.DataFrame({'Group': group})
        temp = pd.concat([group, temp], axis=1)
        fullpeak = fullpeak.append(temp)
    fullpeak = fullpeak.reset_index(drop=True)
    fullpeak.to_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\fullpeak.csv', index=False)
    # Concat successes
    fullsuccess = pd.DataFrame([])
    for date in trials.File.unique():
        # reorganize dataframes
        dframe = pd.read_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\{0}\\{1}_success.csv'.format(
            unicode(date, 'utf-8'), unicode(date, 'utf-8')))
        sdf = pd.concat([dframe['Odor'], dframe[[col for col in dframe.columns if ')_' in col]]], axis=1)
        temp = sdf.T
        temp.reset_index(level=0, inplace=True)
        temp.columns = temp.iloc[0]
        temp = temp.reindex(temp.index.drop(0))
        temp.rename(columns={'Odor': 'Mouse'}, inplace=True)
        temp = temp.reset_index(drop=True)
        group = []
        for x in list(temp.index.values):
            temp.iloc[x]['Mouse'] = temp.iloc[x]['Mouse'].split(')_')[1]
            indexnum = trials.File.unique().tolist().index(temp['Mouse'][x])
            groupname = trials.loc[trials.File == trials.File.unique()[indexnum]].Group.iloc[0]
            group.append(groupname)
        group = pd.DataFrame({'Group': group})
        temp = pd.concat([group, temp], axis=1)
        fullsuccess = fullsuccess.append(temp)
    fullsuccess = fullsuccess.reset_index(drop=True)
    fullsuccess.to_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\fullsuccess.csv', index=False)
    # Concat the integrals
    fullint = pd.DataFrame([])
    for date in trials.File.unique():
        # reorganize dataframes
        dframe = pd.read_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\{0}\\{1}_meanint.csv'.format(
            unicode(date, 'utf-8'), unicode(date, 'utf-8')))
        idf = pd.concat([dframe['Odor'], dframe[[col for col in dframe.columns if ')_' in col]]], axis=1)
        temp = idf.T
        temp.reset_index(level=0, inplace=True)
        temp.columns = temp.iloc[0]
        temp = temp.reindex(temp.index.drop(0))
        temp.rename(columns={'Odor': 'Mouse'}, inplace=True)
        temp = temp.reset_index(drop=True)
        group = []
        for x in list(temp.index.values):
            temp.iloc[x]['Mouse'] = temp.iloc[x]['Mouse'].split(')_')[1]
            indexnum = trials.File.unique().tolist().index(temp['Mouse'][x])
            groupname = trials.loc[trials.File == trials.File.unique()[indexnum]].Group.iloc[0]
            group.append(groupname)
        group = pd.DataFrame({'Group': group})
        temp = pd.concat([group, temp], axis=1)
        fullint = fullint.append(temp)
    fullint = fullint.reset_index(drop=True)
    fullint.to_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\fullintegral.csv', index=False)
    # Get full baseline dataframe
    # Full Baseline
    fullbase = pd.DataFrame([])
    for date in trials.File.unique():
        # reorganize dataframes
        dframe = pd.read_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\{0}\\{1}_meanbase.csv'.format(
            unicode(date, 'utf-8'), unicode(date, 'utf-8')))
        idf = pd.concat([dframe['Odor'], dframe[[col for col in dframe.columns if ')_' in col]]], axis=1)
        temp = idf.T
        temp.reset_index(level=0, inplace=True)
        temp.columns = temp.iloc[0]
        temp = temp.reindex(temp.index.drop(0))
        temp.rename(columns={'Odor': 'Mouse'}, inplace=True)
        temp = temp.reset_index(drop=True)
        group = []
        for x in list(temp.index.values):
            temp.iloc[x]['Mouse'] = temp.iloc[x]['Mouse'].split(')_')[1]
            indexnum = trials.File.unique().tolist().index(temp['Mouse'][x])
            groupname = trials.loc[trials.File == trials.File.unique()[indexnum]].Group.iloc[0]
            group.append(groupname)
        group = pd.DataFrame({'Group': group})
        temp = pd.concat([group, temp], axis=1)
        fullbase = fullbase.append(temp)
    fullbase = fullbase.reset_index(drop=True)
    fullbase.to_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\fullbaseline.csv', index=False)
##