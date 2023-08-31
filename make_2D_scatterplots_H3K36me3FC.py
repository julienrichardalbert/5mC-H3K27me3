#!/usr/bin/env python
# coding: utf-8

# In[1]:


# I'm gonna make a 2D scatterplot with density contour here.
import sys
import seaborn as sb
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None # stop annoying warning when clipping columns. default='warn'


# In[118]:
userInputFile = 'data/10kb_bins_data2.txt'

inputFile = pd.read_csv(userInputFile, sep='\t')
inputFile

x = 'E14_D7_DnmtWT_H3K27me3_CUTnTAG_AMS062021_rep1-2.bam_RPKM'
y = 'E14_D7_DnmtTKO_H3K27me3_CUTnTAG_AMS062021_rep1-2.bam_RPKM'
z = 'K36_logFC_TKO.7_vs_WT.7'

# z = 'E14_D7_DnmtWT_H3K36me3_CnT_AMS052022_rep1-3.bam_RPKM'
#RPKM_MAX = 100          # Determine what is the best max in your dataset
#subset_count = 100000

df = inputFile[[x,y,z]]


# In[144]:


# E14_D0_DnmtWT_H3K27me3_CnT_AMS112021_rep1-2 5.266
# E14_D0_DnmtTKO_H3K27me3_CnT_AMS112021_rep1-2 3.421
# E14_D7_WT_H3K27me3_CUTnTAG_AMS062021_rep1-2 0.865
# E14_D7_TKO_H3K27me3_CUTnTAG_AMS062021_rep1-2 0.298

# normalize by drosophila spike-in ratio
# df[x] = np.log10(df[x] * 0.298 + 0.01)
# df[y] = np.log10(df[y] * 0.895 + 0.01)

# do not normalize by drosophila spike-in ratio
df[x] = np.log10(df[x] + 0.01)
df[y] = np.log10(df[y] + 0.01)


# df = df.dropna(how='any')
df
subset_count = len(df)

# Replace NaN with 0 for the specified columns
df[z] = df[z].fillna(0)


# In[145]:


# OPTIONAL
# Julien now prefers to NOT subsample (set subset_count to number of filtered rows)
# in preference of saving a small but high-res .png file
df_subset = df.sample(n=subset_count)
df_subset


# In[149]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (5,4.1) # makes a perfect square 2D scatterplot (with font size =8 and a 3rd dimension)

fig = df_subset.plot(kind="scatter", x = x, y = y, c = z, \
                     vmin=-2, vmax=2, colormap="RdYlBu_r", alpha=1, s=1, sharex=False)


# In[150]:


sb.kdeplot(x=df[x], y=df[y], cut=0, color="Black", alpha=1)


# In[153]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8

plt.rcParams["figure.figsize"] = (5,4.1) # makes a perfect square 2D scatterplot (with font size=8 and a 3rd dimension)

fig = df_subset.plot(kind="scatter", x = x, y = y, c = z, \
                     vmin=-1, vmax=1, ylim=[-2.15, 2.5], xlim=[-2.15, 2.5], \
                     colormap="RdYlBu_r", alpha=1, s=1, sharex=False)
sb.kdeplot(x=df[x], y=df[y], cut=0, color="Black", alpha=1, levels=5, linewidths=0.8)


# In[154]:


datapoint_count = len(df)
outFigure = "%s_2D_scatter_contour_noNaN_%s_datapoints_%s_subsetShown.svg" % (userInputFile, datapoint_count, subset_count)
fig.figure.savefig(outFigure)
outFigure = "%s_2D_scatter_contour_noNaN_%s_datapoints_%s_subsetShown.png" % (userInputFile, datapoint_count, subset_count)
fig.figure.savefig(outFigure, dpi=600)
