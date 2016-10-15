'''Prep stuff
'C:\Users\Annie\Documents\Data\Ca_Imaging\Analysis\Odor_Panel\Composite_MaxDF_NoP.csv','C:\Users\Annie\Documents\Data\Ca_Imaging\Analysis\Odor_Panel\Composite_Baseline_NoP.csv'


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
from scipy.stats import zscore
from scipy.stats import gaussian_kde
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
def getstd(x):
    '''
    MEAN
    Calculate means for each odor and plots means in boxplot and barplot
    Input: dataframe like comp_sorted
    '''
    Cctrl=x[x['Group']=='Control']
    Cstd=Cctrl.std()
    M=x[x['Group']=='Mint']
    Mstd=M.std()
    H=x[x['Group']=='Hexanal']
    Hstd=H.std()

    CtrlMDFT=pd.DataFrame(Cstd).transpose()
    MMDFT=pd.DataFrame(Mstd).transpose()
    HMDFT=pd.DataFrame(Hstd).transpose()

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
    finalstd=pd.concat(final)
    finalstddf=pd.melt(finalstd,"Group",var_name="Odor")

    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=finalstddf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('Standard Deviation', fontsize=48);
    plt.title('Standard Deviation', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
    return finalstd
##
def getvar(x):
    '''
    Variance (scipy.variation)
    Calculate variation for each odor and plots them, also plots based on concentration
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
    plt.ylabel('Variance', fontsize=48);
    plt.title('Variance Peak DF/F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});

    #calculate and plot variance based on concentration
    MSconcvar = finalvariance[['Group', 'MS 0.01', 'MS 0.05', 'MS 0.1']]
    Hexconcvar = finalvariance[['Group', 'Hexanal 0.01', 'Hexanal 0.05', 'Hexanal 0.1']]
    IAAconcvar = finalvariance[['Group', 'IAA 0.01', 'IAA 0.05', 'IAA 0.1']]

    MSconcvardf = pd.melt(MSconcvar, "Group", var_name="Odor")
    Hexconcvardf = pd.melt(Hexconcvar, 'Group', var_name='Odor')
    IAAconcvardf = pd.melt(IAAconcvar, 'Group', var_name='Odor')

    # Plot variance MS Concentration
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=MSconcvardf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('Variance', fontsize=48);
    plt.title('Variance Peak DF/F - Methyl Salicylate', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});

    # Plot variance Hexanal Concentration
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=Hexconcvardf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('Variance', fontsize=48);
    plt.title('Variance Peak DF/F - Hexanal', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});

    # Plot variance IAA Concentration
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=IAAconcvardf,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('Variance', fontsize=48);
    plt.title('Variance Peak DF/F - Isoamyl acetate', fontsize=55);
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
    #plotcorr.plot(x='DFF', y='Baseline', kind='scatter')
    sns.regplot(plotcorr['DFF'],plotcorr['Baseline']);
    plt.ylabel('Baseline Intensity');
    plt.title('Peak DF/F vs. Baseline');
    plt.xlabel('Peak DF/F');
    from scipy.stats import spearmanr
    rho, pval = spearmanr(plotcorr['DFF'], plotcorr['Baseline'])
    print 'rho, pval'
    print rho, pval
##
def graphAllbox(dataframe):
    #Make Composite or Baseline box plots
    '''
    Argument: DataFrame
    '''
    xmelt=pd.melt(dataframe,'Group',var_name='Odor')
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
def graphAllbar(dataframe):
    '''
    Make Composite or Baseline bargraphs
    Argument: DataFrame
    '''
    xmelt = pd.melt(dataframe, 'Group', var_name='Odor')
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
def concboxplot(odor):
    '''Make box plots based on concentrations'''
    xdf=comp_sorted[['Group','%s01'%odor,'%s05'%odor,'%s10'%odor]]
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
def concbargraph(odor):
    '''Make bar graphs based on concentrations'''
    xdf = comp_sorted[['Group', '%s01' %odor, '%s05' %odor, '%s10' %odor]]
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
def conchist(odor,group):
    '''Make a bunch of plots for an odor ('Mint',etc)'''
    xdf = comp_sorted[['Group', '%s01' % odor, '%s05' % odor, '%s10' % odor]]
    xgroup = xdf[xdf['Group'] == '%s'%group]
    #Set up figure parameters
    sns.set(style="white", palette="muted", color_codes=True)
    sns.set_context("talk", font_scale=2)
    f, axes = plt.subplots(2, 2, figsize=(30, 20), sharex=True)
    f.suptitle('%s group'%group, fontsize=40)
    sns.despine(left=True)
    # data
    d = xgroup['%s01'%odor]
    e = xgroup['%s05'%odor]
    f = xgroup['%s10'%odor]

    # Plot a simple histogram with binsize determined automatically
    sns.distplot(d, kde=False, color="b", hist_kws={"histtype": 'step', "linewidth": 3, "alpha": 0.7}, axlabel=False,
                 ax=axes[0, 0])
    sns.distplot(e, kde=False, color="g", hist_kws={"histtype": 'step', "linewidth": 3, "alpha": 0.7}, axlabel=False,
                 ax=axes[0, 0])
    sns.distplot(f, kde=False, color="r", hist_kws={"histtype": 'step', "linewidth": 3, "alpha": 0.7}, axlabel=False,
                 ax=axes[0, 0])

    # Plot a kernel density estimate and rug plot
    sns.distplot(d, hist=False, rug=True, color="b", axlabel=False, ax=axes[0, 1], label="%s 0.01"%odor)
    sns.distplot(e, hist=False, rug=True, color="g", axlabel=False, ax=axes[0, 1], label="%s 0.05"%odor)
    sns.distplot(f, hist=False, rug=True, color="r", axlabel=False, ax=axes[0, 1], label="%s 0.1"%odor)

    # Plot a filled kernel density estimate
    sns.distplot(d, hist=False, color="b", kde_kws={"shade": True}, axlabel=False, ax=axes[1, 0])
    sns.distplot(e, hist=False, color="g", kde_kws={"shade": True}, axlabel=False, ax=axes[1, 0])
    sns.distplot(f, hist=False, color="r", kde_kws={"shade": True}, axlabel=False, ax=axes[1, 0])

    # Plot a historgram and kernel density estimate
    sns.distplot(d, color="b", axlabel=False, ax=axes[1, 1])
    sns.distplot(e, color="g", axlabel=False, ax=axes[1, 1])
    sns.distplot(f, color="r", axlabel=False, ax=axes[1, 1])

    plt.setp(axes, yticks=[])
    plt.tight_layout()
##
def compare_conc_kruskal(odor):
    '''Do a kruskal wallis test looking at different concentrations of odor
    '''
    xdf = comp_sorted[['Group', '%s01' % odor, '%s05' % odor, '%s10' % odor]]
    xctrl = xdf[xdf['Group'] == 'Control']
    xMS = xdf[xdf['Group'] == 'Mint']
    xHex = xdf[xdf['Group'] == 'Hexanal']
    kctrl=kruskal(xctrl['%s01'%odor],xctrl['%s05'%odor],xctrl['%s10'%odor],nan_policy='omit')
    kmint = kruskal(xMS['%s01' % odor], xMS['%s05' % odor], xMS['%s10' % odor], nan_policy='omit')
    khex = kruskal(xHex['%s01' % odor], xHex['%s05' % odor], xHex['%s10' % odor], nan_policy='omit')
    print 'Control group'
    print kctrl
    print 'Mint group'
    print kmint
    print 'Hexanal group'
    print khex
##
def get_by_odor(odor):
    '''For specified odor, gives boxplot, barplot, histogram, kruskal stat, all across groups
    '''
    x_full=comp_sorted[['Group',odor]]
    xdf=pd.melt(x_full,'Group',var_name='Odor')
    #barplot
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=3);
    plt.figure(figsize=(15, 18));
    sns.barplot(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"}, data=xdf);
    sns.despine();
    plt.legend(loc='upper right');
    plt.ylabel('Peak DF/F');
    plt.title('Peak DF/F');
    #boxplot
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=3);
    plt.figure(figsize=(15, 18));
    sns.boxplot(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"}, data=xdf);
    sns.despine();
    plt.legend(loc='upper right');
    plt.ylabel('Peak DF/F');
    plt.title('Peak DF/F');
    plt.xlabel('Odor');
    #histograms
    xctrl = x_full[x_full['Group'] == 'Control']
    xMS = x_full[x_full['Group'] == 'Mint']
    xHex = x_full[x_full['Group'] == 'Hexanal']
    sns.set(style="white", palette="muted", color_codes=True)
    sns.set_context("talk", font_scale=2)
    # Set up the matplotlib figure
    f, axes = plt.subplots(2, 2, figsize=(30, 20), sharex=True)
    f.suptitle("%s"%odor, fontsize=44)
    sns.despine(left=True)
    # data
    d = xctrl[odor]
    e = xMS[odor]
    f = xHex[odor]
    # Plot a simple histogram with binsize determined automatically
    sns.distplot(d, kde=False, color="r", hist_kws={"histtype": 'step', "linewidth": 3, "alpha": 0.7}, axlabel=False,
                 ax=axes[0, 0])
    sns.distplot(e, kde=False, color="g", hist_kws={"histtype": 'step', "linewidth": 3, "alpha": 0.7}, axlabel=False,
                 ax=axes[0, 0])
    sns.distplot(f, kde=False, color="b", hist_kws={"histtype": 'step', "linewidth": 3, "alpha": 0.7}, axlabel=False,
                 ax=axes[0, 0])
    # Plot a kernel density estimate and rug plot
    sns.distplot(d, hist=False, rug=True, color="r", axlabel=False, ax=axes[0, 1], label="Control")
    sns.distplot(e, hist=False, rug=True, color="g", axlabel=False, ax=axes[0, 1], label="Mint")
    sns.distplot(f, hist=False, rug=True, color="b", axlabel=False, ax=axes[0, 1], label="Hexanal")
    # Plot a filled kernel density estimate
    sns.distplot(d, hist=False, color="r", kde_kws={"shade": True}, axlabel=False, ax=axes[1, 0])
    sns.distplot(e, hist=False, color="g", kde_kws={"shade": True}, axlabel=False, ax=axes[1, 0])
    sns.distplot(f, hist=False, color="b", kde_kws={"shade": True}, axlabel=False, ax=axes[1, 0])
    # Plot a historgram and kernel density estimate
    sns.distplot(d, color="r", axlabel=False, ax=axes[1, 1])
    sns.distplot(e, color="g", axlabel=False, ax=axes[1, 1])
    sns.distplot(f, color="b", axlabel=False, ax=axes[1, 1])
    plt.setp(axes, yticks=[])
    plt.tight_layout()
    # Stats
    print kruskal(xctrl['%s' % odor], xMS['%s' % odor], xHex['%s' % odor], nan_policy='omit')
##
def baseline_sort_by_group(baselinedf):
    '''Sort baseline by group and make some graphs
    use variable bdf from getdata function'''
    bdfull = pd.melt(baselinedf, "Group", var_name="Odor")
    bdfull = bdfull.dropna()
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=2.2);
    plt.figure(figsize=(45, 20));
    ax = sns.barplot(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"},
                     data=bdfull);
    sns.despine()
    plt.ylabel('Intensity', fontsize=48);
    plt.title('Baseline F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 30});
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=2.2);
    plt.figure(figsize=(45, 20));
    ax = sns.boxplot(x="Odor", y="value", hue="Group", palette={"Control": "r", "Hexanal": "b", "Mint": "g"},
                     data=bdfull);
    sns.despine()
    plt.ylabel('Intensity', fontsize=48);
    plt.title('Baseline F', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 30});
##
def rank_odor(dataframe):
    '''Rank each odor for each individual cell'''
    x=dataframe.rank(axis=1,numeric_only=float,na_option='keep',ascending=False)
    x_label=pd.DataFrame(dataframe.Group)
    tmp=[x_label,x]
    rankedx=pd.concat(tmp,axis=1)
    ctrlmr=pd.DataFrame(rankedx[rankedx['Group']=='Control'].mean(axis=0)).T
    mintmr=pd.DataFrame(rankedx[rankedx['Group']=='Mint'].mean(axis=0)).T
    hexmr=pd.DataFrame(rankedx[rankedx['Group']=='Hexanal'].mean(axis=0)).T
    all=pd.concat([ctrlmr,mintmr,hexmr])
    all=all.reset_index(drop=True)
    grouplabels=pd.DataFrame({'Group':['Control','Mint','Hexanal']})
    finaldf=pd.concat([grouplabels,all],axis=1)
    finalmelt = pd.melt(finaldf, "Group", var_name="Odor")
    #Plot some stuff
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(55, 20));
    sns.pointplot(x="Odor", y="value", hue="Group", data=finalmelt,
                  palette={"Control": "r", "Mint": "g", 'Hexanal': 'b'});
    sns.despine()
    plt.ylabel('Rank', fontsize=48);
    plt.title('Rank of odor response', fontsize=55);
    plt.xlabel('Odor', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
    print finaldf
    print kruskal(ctrlmr,mintmr,hexmr)
##
def zodors(dataframe,odor):
    '''Compare the zscores of the odor dataframes, based on odor and comparison across groups'''
    df=dataframe[['Group',odor]]
    ctrl=zscore(df[df['Group']=='Control'][odor])
    mint=zscore(df[df['Group']=='Mint'][odor])
    hex=zscore(df[df['Group']=='Hexanal'][odor])
    #plot some stuff
    sns.set(style="white", palette="muted", color_codes=True)
    plt.figure(figsize=(24, 18));
    sns.distplot(ctrl, color="r", hist=False,kde_kws={"shade": True},axlabel=False, label="Control")
    sns.distplot(mint, color="g", hist=False,kde_kws={"shade": True},axlabel=False, label="Mint")
    sns.distplot(hex, color="b", hist=False,kde_kws={"shade": True},axlabel=False, label="Hexanal")
    sns.despine();
    plt.legend(loc='upper right');
    plt.xlabel('Z-score');
    plt.title('Z-scores for %s'%odor);
    return anderson_ksamp((ctrl,mint,hex))
##
def kdeshift(x,odor):
    '''Compare kde visually by shifting all kdes
    so that all peaks are at zero
    x is the dataframe (generally peaks.comp_sorted)
    '''
    ctrl = x[x['Group'] == 'Control']
    ms = x[x['Group'] == 'Mint']
    hex = x[x['Group'] == 'Hexanal']
    x=np.arange(-1,5,0.01)
    # Get the peak of each kde
    ckde = gaussian_kde(ctrl[odor])
    cy = ckde.evaluate(x)
    mkde = gaussian_kde(ms[odor])
    my = mkde.evaluate(x)
    hkde = gaussian_kde(hex[odor])
    hy = hkde.evaluate(x)
    # use this max x value to shift the entire distribution
    c_max_x_value = x[np.argmax(cy)]
    m_max_x_value = x[np.argmax(my)]
    h_max_x_value = x[np.argmax(hy)]
    # make a dataframe of x,y coordinates of the kde, with a bin size of 0.01
    cxdf = pd.DataFrame(x)
    cxdf.columns = ['x']
    c_newxdf = pd.DataFrame(cxdf['x'] - c_max_x_value)
    cydf = pd.DataFrame(cy)
    cydf.columns = ['y']
    c_coordinates = pd.concat([cxdf, cydf], axis=1)
    c_newcoordinates = pd.concat([c_newxdf, cydf], axis=1)

    mxdf = pd.DataFrame(x)
    mxdf.columns = ['x']
    m_newxdf = pd.DataFrame(mxdf['x'] - m_max_x_value)
    mydf = pd.DataFrame(my)
    mydf.columns = ['y']
    m_coordinates = pd.concat([mxdf, mydf], axis=1)
    m_newcoordinates = pd.concat([m_newxdf, mydf], axis=1)

    hxdf = pd.DataFrame(x)
    hxdf.columns = ['x']
    h_newxdf = pd.DataFrame(hxdf['x'] - h_max_x_value)
    hydf = pd.DataFrame(hy)
    hydf.columns = ['y']
    h_coordinates = pd.concat([hxdf, hydf], axis=1)
    h_newcoordinates = pd.concat([h_newxdf, hydf], axis=1)

    # plot it!
    sns.set(palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=3);
    plt.figure(figsize=(24, 18));
    plt.plot(c_newcoordinates['x'], c_newcoordinates['y'], color='r', label='Control');
    plt.plot(m_newcoordinates['x'], m_newcoordinates['y'], color='g', label='Mint');
    plt.plot(h_newcoordinates['x'], h_newcoordinates['y'], color='b', label='Hexanal');
    sns.despine();
    plt.legend(loc='upper right');
    plt.title('KDEs, peaks centered, %s'%odor);
##
def hist_analysis(dataframe,odor):
    '''Compare normalized histograms'''
    xfull = dataframe[['Group', odor]]
    xdf = pd.melt(xfull, "Group", var_name="Odor")
    cg = xfull[xfull['Group'] == 'Control']
    mg = xfull[xfull['Group'] == 'Mint']
    hg = xfull[xfull['Group'] == 'Hexanal']
    #Make histograms
    ctrlhist = np.asarray(np.histogram(cg[odor], bins=140, range=(-1, 6))[0], dtype=np.float)
    mshist = np.asarray(np.histogram(mg[odor], bins=140, range=(-1, 6))[0], dtype=np.float)
    hexhist = np.asarray(np.histogram(hg[odor], bins=140, range=(-1, 6))[0], dtype=np.float)
    ctrlnorm = ctrlhist / (ctrlhist.sum())
    msnorm = mshist / (mshist.sum())
    hexnorm = hexhist / (hexhist.sum())
    bincounts = np.histogram(cg[odor], bins=140, range=(-1, 6))[1][0:-1]
    # dataframe of bin values
    colnames = ['Bin', 'Control', 'Mint', 'Hexanal']
    tmp = [pd.DataFrame(bincounts), pd.DataFrame(ctrlnorm), pd.DataFrame(msnorm), pd.DataFrame(hexnorm)]
    fullnorm = pd.concat(tmp, axis=1)
    fullnorm.columns = colnames
    # plot the bins
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(30, 15));
    plt.scatter(fullnorm['Bin'], fullnorm['Control'], s=80, c='r', label='Control');
    plt.scatter(fullnorm['Bin'], fullnorm['Mint'], s=80, c='g', label='Mint');
    plt.scatter(fullnorm['Bin'], fullnorm['Hexanal'], s=80, c='b', label='Hexanal');
    plt.xlim(-0.1, 4.1)
    sns.despine()
    plt.ylabel('Proportion', fontsize=48);
    plt.title('Normalized Histograms', fontsize=55);
    plt.xlabel('Peak DF/F', fontsize=48);
    plt.legend(loc=1, prop={'size': 48});
    #Get some relevant values
    print 'Kruskal-Wallis of %s between exposure groups'%odor
    print kruskal(fullnorm['Control'],fullnorm['Mint'],fullnorm['Hexanal'])
    print 'Kruskal-Wallis of Hexanal and Control'
    print kruskal(fullnorm['Control'],fullnorm['Hexanal'])
    print 'KS test of Hexanal and Control'
    print ks_2samp(fullnorm['Hexanal'],fullnorm['Control'])
    print 'Kruskal-Wallis of Hexanal and Mint'
    print kruskal(fullnorm['Hexanal'],fullnorm['Mint'])
    print 'KS test of Hexanal and Mint'
    print ks_2samp(fullnorm['Hexanal'],fullnorm['Mint'])
    print 'Kruskal-Wallis of Mint and Control'
    print kruskal(fullnorm['Mint'],fullnorm['Control'])
    print 'KS test of Mint and Control'
    print ks_2samp(fullnorm['Mint'], fullnorm['Control'])
##
def nsfa_by_group(dataframe):
    '''Non stationary fluctuation analysis, averaged by group,
    not averaged between odors
    '''
    # Calculate means and variance for each odor
    Cctrl = dataframe[dataframe['Group'] == 'Control']
    Cmean = pd.DataFrame(Cctrl.mean())
    Cmean.columns = ['Control Mean']
    Cvar = pd.DataFrame(Cctrl.var())
    Cvar.columns = ['Control Variance']
    M = dataframe[dataframe['Group'] == 'Mint']
    Mmean = pd.DataFrame(M.mean())
    Mmean.columns = ['Mint Mean']
    Mvar = pd.DataFrame(M.var())
    Mvar.columns = ['Mint Variance']
    H = dataframe[dataframe['Group'] == 'Hexanal']
    Hmean = pd.DataFrame(H.mean())
    Hmean.columns = ['Hexanal Mean']
    Hvar = pd.DataFrame(H.var())
    Hvar.columns = ['Hexanal Variance']
    # Concat
    Ctmp = [Cmean, Cvar]
    Mtmp = [Mmean, Mvar]
    Htmp = [Hmean, Hvar]
    CtrlDF = pd.concat(Ctmp, axis=1)
    MDF = pd.concat(Mtmp, axis=1)
    HDF = pd.concat(Htmp, axis=1)
    final = [CtrlDF, MDF, HDF]
    finaldf = pd.concat(final, axis=1)
    finaldf = finaldf.reset_index(drop=True)
    #Plot it
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(30, 15));
    sns.regplot(finaldf['Control Mean'], finaldf['Control Variance'], scatter_kws={"s": 175}, color='r',label='Control')
    sns.regplot(finaldf['Mint Mean'], finaldf['Mint Variance'], scatter_kws={"s": 175}, color='g',label='Mint')
    sns.regplot(finaldf['Hexanal Mean'], finaldf['Hexanal Variance'], scatter_kws={"s": 175}, color='b',label='Hexanal')
    sns.despine()
    plt.ylabel('Variance', fontsize=48);
    plt.title('Mean vs. Variance', fontsize=55);
    plt.xlabel('Mean', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
##
def nsfa_by_cell(dataframe):
    '''Non stationary fluctuation analysis, averaged by cell across odors'''
    Cctrl = dataframe[dataframe['Group'] == 'Control']
    M = dataframe[dataframe['Group'] == 'Mint']
    H = dataframe[dataframe['Group'] == 'Hexanal']
    Ccellmean = Cctrl.mean(axis=1)
    Ccellvar = Cctrl.var(axis=1)
    Mcellmean = M.mean(axis=1)
    Mcellvar = M.var(axis=1)
    Hcellmean = H.mean(axis=1)
    Hcellvar = H.var(axis=1)
    # Concat
    Ctemp = [Cctrl['Group'], Ccellmean, Ccellvar]
    Mtemp = [M['Group'], Mcellmean, Mcellvar]
    Htemp = [H['Group'], Hcellmean, Hcellvar]
    CtrlcellDF = pd.concat(Ctemp, axis=1)
    CtrlcellDF.columns = ('Group', 'Mean', 'Variance')
    McellDF = pd.concat(Mtemp, axis=1)
    McellDF.columns = ('Group', 'Mean', 'Variance')
    HcellDF = pd.concat(Htemp, axis=1)
    HcellDF.columns = ('Group', 'Mean', 'Variance')
    finalcell = [CtrlcellDF, McellDF, HcellDF]
    finalcelldf = pd.concat(finalcell, axis=0)
    sns.set(style="white", palette="muted", color_codes=True);
    sns.set_context("talk", font_scale=1.8);
    plt.figure(figsize=(30, 15));
    sns.regplot('Mean', 'Variance', CtrlcellDF, scatter_kws={"s": 80}, color='r', label='Control')
    sns.regplot('Mean', 'Variance', McellDF, scatter_kws={"s": 80}, color='g', label='Mint')
    sns.regplot('Mean', 'Variance', HcellDF, scatter_kws={"s": 80}, color='b', label='Hexanal')
    sns.despine()
    plt.ylabel('Variance', fontsize=48);
    plt.title('Mean vs. Variance', fontsize=55);
    plt.xlabel('Mean', fontsize=48);
    plt.legend(loc=2, prop={'size': 48});
##