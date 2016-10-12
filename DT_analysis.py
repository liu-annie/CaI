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
def getpeaks(session):
    '''Spits out all the peaks from imaging session
    session input as string
    '''
    temp=pd.DataFrame([])
    os.chdir('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\%s'%session)
    for filename in glob.glob('*dt.txt'):
        df=pd.read_csv(filename,nrows=175)
        df=df[[col for col in df.columns if 'G PMT' in col]]
        peak=[]
        for col in df.columns:
            a=df[col]
            firsta=1;
            firstb=24;
            #Figures out if there is a min or max and sees if it passes threshold (3SD)
            if np.absolute(min(a[26:80]))>np.absolute(max(a[26:80])) and np.absolute(min(a[26:80]))>=3*np.std(df[col][firsta:firstb]):
                b=min(a[26:80])
                peak.append(b)
            elif np.absolute(max(a[26:80]))>np.absolute(min(a[26:80]))and np.absolute(max(a[26:80]))>=3*np.std(df[col][firsta:firstb]):
                b=max(a[26:80])
                peak.append(b)
            else:
                b=0
                peak.append(b)
        peaks=pd.DataFrame(peak).T
        temp=temp.append(peaks,ignore_index=True)
    temp.columns=df.columns
    peakdf=pd.concat([pd.DataFrame({'Trial':np.arange(1,len(temp)+1,1)}),temp],axis=1)
    peakdf.to_csv('%s_peakdf.txt' % session)
    #Get the odor trial order
    trials = pd.read_csv('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\Analysis\\Odor_Panel\\Odor_Trials.csv')
    filerow = trials.loc[trials['File'] == session]
    odortrials = {}
    for t in trials.Odor.unique():
        y = {t: [int(x) for x in filerow.loc[filerow['Odor'] == t][['T1', 'T2', 'T3', 'T4']].values.tolist()[0]]}
        odortrials.update(y)
    # Get average peak across all trials using peakdf dataframe
    meandf = pd.DataFrame([])
    for key in odortrials:
        odor = odortrials[key]
        peakdf.loc[peakdf['Trial'].isin(odor)]
        mean = []
        for col in peakdf.loc[peakdf['Trial'].isin(odor)][
            [col for col in peakdf.loc[peakdf['Trial'].isin(odor)].columns if 'G PMT' in col]]:
            mean.append(peakdf.loc[peakdf['Trial'].isin(odor)][col].mean())
        mean = pd.DataFrame(mean).T
        mean.columns = peakdf.loc[peakdf['Trial'].isin(odor)][
            [col for col in peakdf.loc[peakdf['Trial'].isin(odor)].columns if 'G PMT' in col]].columns
        meandf = meandf.append(mean)
    meandf = meandf.reset_index(drop=True)
    meandf.columns = [str(col) + '_' + session for col in meandf.columns]
    meandf = pd.concat([pd.DataFrame({'Odor': odortrials.keys()}), meandf], axis=1)
    meandf.to_csv('%s_mean.txt' % session)
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
    successdf.columns = [str(col) + '_' + session for col in successdf.columns]
    successdf = pd.concat([pd.DataFrame({'Odor': odortrials.keys()}), successdf], axis=1)
    successdf.to_csv('%s_success.txt' % session)