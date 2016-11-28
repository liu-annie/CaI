#Dependencies
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
def detrend(file):
    '''script to detrend files and spit out notepad files of detrended data
    'C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\GoodFiles\\160321_2'
    '''
    os.chdir('C:\\Users\\Annie\\Documents\\Data\\Ca_Imaging\\HabituationFiles\\%s'%file)
    for filename in glob.glob('*.txt'):
        df=pd.read_table(filename,skiprows=4)
        odf=df[[col for col in df.columns if 'G PMT' in col]]
        # dff is DF/F for every cell
        base=[];
        for col in odf.columns:
            value=odf.loc[2:23,col].mean();
            base.append(value)
        base=pd.DataFrame(base).T
        base.columns=odf.columns
        # base
        dff=[];
        for col in odf.columns:
            value=(odf[col]-np.asscalar(base[col]))/np.asscalar(base[col])
            dff.append(value)
        dff=pd.DataFrame(dff).T
        # detrend
        dt=[];
        firsta=1;
        firstb=20;
        lasta=130
        lastb=174
        for col in dff.columns:
            firstvalues=pd.DataFrame([np.arange(firsta,firstb+1,1),dff.loc[firsta:firstb,col]]).T
            lastvalues=pd.DataFrame([np.arange(lasta,lastb+1,1),dff.loc[lasta:lastb,col]]).T
            temp=firstvalues.append(lastvalues)
            fitline=np.polyfit(temp[0],temp[1],2)
            tmp=np.poly1d(fitline)
            fitvalues=tmp(np.arange(len(dff[col])))
            detrended=np.subtract(dff[col],fitvalues)
            dt.append(detrended)
        dt=pd.DataFrame(dt).T
        dt.to_csv(filename.split('.')[0]+'dt.txt')