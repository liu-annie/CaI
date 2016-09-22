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
%matplotlib
import seaborn as sns
##
'''Import and organize data'''
##
comp=pd.read_csv('C:\Users\Annie\Documents\Data\Ca_Imaging\Analysis\Odor_Panel\Composite_MaxDF_NoP.csv')
del comp['Mouse']
comp_sorted=comp.reindex_axis(comp.mean().sort_values().index, axis=1)
comp_labels=pd.DataFrame(comp.Group)
tmp=[comp_labels,comp_sorted]
composite_full=pd.concat(tmp,axis=1)
composite_full.head()
cfull=pd.melt(composite_full,"Group",var_name="Odor")
##
'''Composite Analyses'''
##
#Calculate means for each odor
Cctrl=composite_full[composite_full['Group']=='Control']
Cmean=Cctrl.mean()
M=composite_full[composite_full['Group']=='Mint']
Mmean=M.mean()
H=composite_full[composite_full['Group']=='Hexanal']
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
##

