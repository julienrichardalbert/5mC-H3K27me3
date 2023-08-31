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


# In[2]:
# 10kb bins
userInputFile = 'data/10kb_bins_data.txt'
# CGIs
# userInputFile = 'data/CGI_data.txt'
inputFile = pd.read_csv(userInputFile, sep='\t')
inputFile

# get column names to plot
for column_name in inputFile.columns:
    print(column_name)


# In[19]:
# E14_D0_DnmtWT_H3K27me3_CnT_AMS112021_rep1-2.bam_RPKM
# E14_D0_DnmtTKO_H3K27me3_CnT_AMS112021_rep1-2.bam_RPKM
# E14_D7_DnmtWT_H3K27me3_CUTnTAG_AMS062021_rep1-2.bam_RPKM
# E14_D7_DnmtTKO_H3K27me3_CUTnTAG_AMS062021_rep1-2.bam_RPKM
# E14_D0_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5_mean
# E14_D7_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5_mean

# change which column you want to plot by setting the x,y and z variables here
x = 'E14_D0_DnmtTKO_H3K27me3_CnT_AMS112021_rep1-2.bam_RPKM'
y = 'E14_D7_DnmtTKO_EzhipKO_H3K27me3_CnT_202212_rep1-2.bam_RPKM'
z = 'E14_D7_WGBS_SEPE_raw_rmDup.CpG_report_mergeTwoStrands_x5_mean'
RPKM_MAX = 100          # Determine what is the best max in your dataset
subset_count = 100000

df = inputFile[[x,y,z]]
df


# In[20]:


# do not normalize using drosophila spike-in. headache, similar results
df[x] = np.log10(df[x] + 0.01)
df[y] = np.log10(df[y] + 0.01)

df = df.dropna(how='any')
df
# do not subset
subset_count = len(df)


# In[21]:


# OPTIONAL
# Julien now prefers to NOT subsample (set subset_count to number of filtered rows)
# in preference of saving a small but high-res .png file
df_subset = df.sample(n=subset_count)
df_subset


# In[22]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (5,4.1) # makes a perfect square 2D scatterplot (with font size =8 and a 3rd dimension)

fig = df_subset.plot(kind="scatter", x = x, y = y, c = z, \
                     vmin=0, vmax=100, colormap="RdYlBu_r", alpha=1, s=1, sharex=False)

#fig.set_yscale('log')
#fig.set_xscale('log')


# In[11]:


sb.kdeplot(x=df[x], y=df[y], cut=0, color="Black", alpha=1)


# In[23]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8

plt.rcParams["figure.figsize"] = (5,4.1) # makes a perfect square 2D scatterplot (with font size=8 and a 3rd dimension)

fig = df_subset.plot(kind="scatter", x = x, y = y, c = z, \
                     vmin=0, vmax=100, ylim=[-2.15, 2.5], xlim=[-2.15, 2.5], \
                     colormap="RdYlBu_r", alpha=1, s=1, sharex=False)
sb.kdeplot(x=df[x], y=df[y], cut=0, color="Black", alpha=1, levels=5, linewidths=0.8)


# In[28]:


datapoint_count = len(df)
outFigure = "%s_2D_scatter_contour_noNaN_%s_RPKMmax_%s_datapoints_%s_subsetShown.svg" % (userInputFile, RPKM_MAX, datapoint_count, subset_count)
fig.figure.savefig(outFigure)
outFigure = "%s_2D_scatter_contour_noNaN_%s_RPKMmax_%s_datapoints_%s_subsetShown.png" % (userInputFile, RPKM_MAX, datapoint_count, subset_count)
fig.figure.savefig(outFigure, dpi=600)
