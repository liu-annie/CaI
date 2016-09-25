'''Prep stuff
'''
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
from scipy.stats import spearmanr
import seaborn as sns
##
def getdata(dff_file,baseline_file,exclude):
    '''Import data
    reads csv file, excludes uncessary columns, sorts based on column mean, ascending
    returns 2 global variables that can be used in all subsequent functions
    comp_sorted: composite dataframe sorted based on column mean
    bdf: baseline data frame
    '''
    file=pd.read_csv(dff_file)
    cols=[col for col in file.columns if col not in exclude]
    comp=file[cols]
    comp2=comp.reindex_axis(comp.mean().sort_values().index,axis=1)
    comp_labels=pd.DataFrame(comp.Group)
    tmp=[comp_labels,comp2]
    global comp_sorted
    comp_sorted=pd.concat(tmp,axis=1)

    blfile = pd.read_csv(baseline_file)
    cols = [col for col in blfile.columns if col not in exclude]
    bl = blfile.reindex_axis(blfile.mean().sort_values().index, axis=1)
    bl_labels = pd.DataFrame(blfile.Group)
    tmp = [bl_labels, bl]
    global bdf
    bdf = pd.concat(tmp, axis=1)

##
'''Composite Analyses'''
def getmean(x):
    '''
    MEAN
    Calculate means for each odor and plots means in boxplot and barplot
    Input: dataframe like comp_sorted
    '''
    Cctrl=x[x['Group']=='Control']
    Cmean=Cctrl.mean()
    M=x[x['Group']=='Mint']
    Mmean=M.mean()
    H=x[x['Group']=='Hexanal']
    Hmean=H.mean()

    CtrlMDFT=pd.DataFrame(Cmean).transpose()
    MMDFT=pd.DataFrame(Mmean).transpose()
    HMDFT=pd.DataFrame(Hmean).transpose()

    #add group labels back
    gnc = pd.DataFrame({'Group':['Control']})
    gnm=pd.DataFrame({'Group':['Mint']})
    gnh=pd.DataFrame({'Group':['Hexanal']})

    Ctmp=[gnc,CtrlMDFT]
    Mtmp=[gnm,MMDFT]
    Htmp=[gnh,HMDFT]

    CtrlMDF=pd.concat(Ctmp,axis=1)
    MMDF=pd.concat(Mtmp,axis=1)
    HMDF=pd.concat(Htmp,axis=1)

    final=[CtrlMDF,MMDF,HMDF]
    finalmean=pd.concat(final)
    finalmeandf=pd.melt(finalmean,"Group",var_name="Odor")

    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=finalmeandf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('DF/F', fontsize=48);
    plt.title('Mean Peak DF/F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
    return finalmean
##
def getmedian(x):
    '''
    MEDIAN
    Calculate medians for each odor and plots means in boxplot and barplot
    Input: dataframe like comp_sorted
    '''
    # Calculate medians for each odor
    Cctrl = x[x['Group'] == 'Control']
    Cmedian = Cctrl.median()
    M = x[x['Group'] == 'Mint']
    Mmedian = M.median()
    H = x[x['Group'] == 'Hexanal']
    Hmedian = H.median()

    CtrlMedT = pd.DataFrame(Cmedian).transpose()
    MMedT = pd.DataFrame(Mmedian).transpose()
    HMedT = pd.DataFrame(Hmedian).transpose()

    # add group labels back
    gnc = pd.DataFrame({'Group': ['Control']})
    gnm = pd.DataFrame({'Group': ['Mint']})
    gnh = pd.DataFrame({'Group': ['Hexanal']})

    Ctmp = [gnc, CtrlMedT]
    Mtmp = [gnm, MMedT]
    Htmp = [gnh, HMedT]

    CtrlMed = pd.concat(Ctmp, axis=1)
    MMed = pd.concat(Mtmp, axis=1)
    HMed = pd.concat(Htmp, axis=1)

    finalmed = [CtrlMed, MMed, HMed]
    finalmedian = pd.concat(finalmed)
    finalmediandf = pd.melt(finalmedian, "Group", var_name="Odor")
    # Plot medians
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=finalmediandf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('DF/F', fontsize=48);
    plt.title('Median Peak DF/F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
    return finalmedian
##
def getvar(x):
    '''
    COV
    Calculate COV for each odor and plots them, also plots based on concentration
    Input: dataframe like comp_sorted
    '''
    # Calculate variance
    CC = x[x['Group'] == 'Control']
    MM = x[x['Group'] == 'Mint']
    HH = x[x['Group'] == 'Hexanal']
    del CC['Group']
    del MM['Group']
    del HH['Group']

    CvarT = pd.DataFrame(variation(CC)).transpose()
    MvarT = pd.DataFrame(variation(MM)).transpose()
    HvarT = pd.DataFrame(variation(HH)).transpose()

    # add group labels back
    gnc = pd.DataFrame({'Group': ['Control']})
    gnm = pd.DataFrame({'Group': ['Mint']})
    gnh = pd.DataFrame({'Group': ['Hexanal']})

    Ctmp = [gnc, CvarT]
    Mtmp = [gnm, MvarT]
    Htmp = [gnh, HvarT]

    CV = pd.concat(Ctmp, axis=1)
    MV = pd.concat(Mtmp, axis=1)
    HV = pd.concat(Htmp, axis=1)

    gn = ['Group']
    gname = np.array(gn)
    odornames = CC.columns.values
    columnnames = np.append(gname, odornames)
    CV.columns = columnnames
    MV.columns = columnnames
    HV.columns = columnnames

    # Fix BLANK NaNs
    CCNB = pd.DataFrame(CC['BLANK'])
    CCNB = CCNB.dropna()
    CCNBvar = variation(CCNB)
    CCNBvars = np.asscalar(CCNBvar)
    CV = CV.replace(CV.BLANK[0], CCNBvars)

    MMNB = pd.DataFrame(MM['BLANK'])
    MMNB = MMNB.dropna()
    MMNBvar = variation(MMNB)
    MMNBvars = np.asscalar(MMNBvar)
    MV = MV.replace(MV.BLANK[0], MMNBvars)

    HHNB = pd.DataFrame(HH['BLANK'])
    HHNB = HHNB.dropna()
    HHNBvar = variation(HHNB)
    HHNBvars = np.asscalar(HHNBvar)
    HV = HV.replace(HV.BLANK[0], HHNBvars)

    complete = [CV, HV, MV]
    global finalvariance
    finalvariance = pd.concat(complete)
    finalvariancedf = pd.melt(finalvariance, "Group", var_name="Odor")

    #plot variance
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=finalvariancedf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('COV', fontsize=48);
    plt.title('COV Peak DF/F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});

    #calculate and plot variance based on concentration
    MSconcvar = finalvariance[['Group', 'MS 0.01', 'MS 0.05', 'MS 0.1']]
    Hexconcvar = finalvariance[['Group', 'Hexanal 0.01', 'Hexanal 0.05', 'Hexanal 0.1']]
    IAAconcvar = finalvariance[['Group', 'IAA 0.01', 'IAA 0.05', 'IAA 0.1']]

    MSconcvardf = pd.melt(MSconcvar, "Group", var_name="Odor")
    Hexconcvardf = pd.melt(Hexconcvar, 'Group', var_name='Odor')
    IAAconcvardf = pd.melt(IAAconcvar, 'Group', var_name='Odor')

    # Plot COV MS Concentration
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=MSconcvardf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('COV', fontsize=48);
    plt.title('COV Peak DF/F - Methyl Salicylate', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});

    # Plot COV Hexanal Concentration
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=Hexconcvardf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('COV', fontsize=48);
    plt.title('COV Peak DF/F - Hexanal', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});

    # Plot COV IAA Concentration
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=IAAconcvardf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('COV', fontsize=48);
    plt.title('COV Peak DF/F - Isoamyl acetate', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
##
def corrBLDFF(baseline,dff):
    '''Correlate baseline with DFF
    '''
    #make sure all columns are in the same order
    dff_cn = dff.columns.tolist()
    bldf = baseline[dff_cn]
    # Turn both dataframes into single columns
    bdfmelt = pd.melt(bldf, "Group", var_name='Odor')
    cfmelt = pd.melt(dff, 'Group', var_name='Odor')
    # Make a dataframe so you can plot the correlation
    plotcorr = pd.concat([cfmelt, bdfmelt['value']], axis=1)
    # Rename these columns to what they're supposed to be
    plotcorr.columns = ['Group', 'Odor', 'DFF', 'Baseline']
    # Plot Baseline vs. DFF as a scatterplot
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.4);
    plotcorr.plot(x='DFF', y='Baseline', kind='scatter')
    plt.ylabel('Baseline Intensity');
    plt.title('Peak DF/F vs. Baseline');
    plt.xlabel('Peak DF/F');
    from scipy.stats import spearmanr
    rho, pval = spearmanr(plotcorr['DFF'], plotcorr['Baseline'])
    return rho, pval
##
def graphAllbox(x):
    '''
    Make Composite or Baseline box plots
    Argument: DataFrame
    '''
    xmelt=pd.melt(x,'Group',var_name='Odor')
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=2.2);
    plt.figure(figsize=(45, 20));
    ax = sns.barplot(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"},
                     data=xmelt);
    sns.despine()
    plt.ylabel('Peak DF/F', fontsize=48);
    plt.title('Peak DF/F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
##
def graphAllbar(x):
    '''
    Make Composite or Baseline bargraphs
    Argument: DataFrame
    '''
    xmelt = pd.melt(x, 'Group', var_name='Odor')
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=2.2);
    plt.figure(figsize=(45, 20));
    ax = sns.boxplot(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"},
                     data=xmelt);
    sns.despine()
    plt.ylabel('Peak DF/F', fontsize=48);
    plt.title('Peak DF/F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
##
def concboxplot(x):
    '''Make box plots based on concentrations'''
    xdf=comp_sorted[['Group','%s 0.01'%x,'%s 0.05'%x,'%s 0.1'%x]]
    xdf=pd.melt(xdf,'Group',var_name='Odor')
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=4.3);
    plt.figure(figsize=(30, 20));
    sns.boxplot(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"}, data=xdf);
    sns.despine()
    plt.legend(loc='upper right');
    plt.ylabel('Peak DF/F');
    plt.title('Peak DF/F');
##
def concbargraph(x):
    '''Make bar graphs based on concentrations'''
    xdf = comp_sorted[['Group', '%s 0.01' % x, '%s 0.05' % x, '%s 0.1' % x]]
    xdf = pd.melt(xdf, 'Group', var_name='Odor')
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=4.3);
    plt.figure(figsize=(30, 20));
    sns.bargraph(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"}, data=xdf);
    sns.despine()
    plt.legend(loc='upper right');
    plt.ylabel('Peak DF/F');
    plt.title('Peak DF/F');
##
def conchist(x):
    xctrl = MS_full[MS_full['Group'] == 'Control']
    xMS = MS_full[MS_full['Group'] == 'Mint']
    xhex = MS_full[MS_full['Group'] == 'Hexanal']
    xctrl.tail()