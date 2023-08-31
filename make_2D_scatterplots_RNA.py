#!/usr/bin/env python
# coding: utf-8

# In[2]:


# I'm gonna make a 2D scatterplot with density contour here.
import sys
import seaborn as sb
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
pd.options.mode.chained_assignment = None # stop annoying warning when clipping columns. default='warn'


# In[3]:
userInputFile = 'data/mm10_refseq_genes_and_isoforms_scatterplot.txt'

inputFile = pd.read_csv(userInputFile, sep='\t')
inputFile


WT  = ['E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10.bam_RPKM',  'E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10.bam_RPKM']
TKO = ['E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10.bam_RPKM']


#WT = ['E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10.bam_RPKM']
#TKO = ['E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10.bam_RPKM']

df = inputFile[WT + TKO]
#RPKM_MAX = 100          # Determine what is the best max in your dataset
#subset_count = 100000
df['log_WT'] = np.log10(df[WT].mean(axis=1)+0.01)
df['log_TKO'] = np.log10(df[TKO].mean(axis=1)+0.01)
df['z'] = inputFile['index2']




def color_mapping(values):
    unique_values = values.unique()
    num_unique = len(unique_values)
    cmap = cm.get_cmap('viridis_r', num_unique)  # Choose a colormap based on the number of unique values

    color_dict = {}
    for i, value in enumerate(unique_values):
        color_dict[value] = cmap(i)

    return color_dict

# Example usage:
df['w'] = df['z'].map(color_mapping(df['z']))


# In[5]:


#df = df.dropna(how='any')
df
subset_count = len(df)


# In[6]:


# OPTIONAL
# Julien now prefers to NOT subsample (set subset_count to number of filtered rows)
# in preference of saving a small but high-res .png file
df_subset = df.sample(n=subset_count)
df_subset


# In[7]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (5,4.1) # makes a perfect square 2D scatterplot (with font size =8 and a 3rd dimension)

fig = df_subset.plot(kind="scatter", x = 'log_WT', y = 'log_TKO', c = 'w', \
                     vmin=0, vmax=1, alpha=1, s=25, sharex=False)





# In[8]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (4.1,4.1) # makes a perfect square 2D scatterplot (with font size=8 and a 3rd dimension)

fig = df_subset.plot(kind="scatter", x = 'log_WT', y = 'log_TKO', c = 'w', \
                     vmin=0, vmax=1, ylim=[-2.15, 3.5], xlim=[-2.15, 3.5], \
                     alpha=1, s=22, sharex=False)

sb.kdeplot(x=df['log_WT'], y=df['log_TKO'], cut=0, color="Black", alpha=1, levels=5, linewidths=0.8)


# In[9]:


datapoint_count = len(df)
outFigure = "%s_2D_scatter_contour_%s_datapoints_%s_subsetShown.svg" % (userInputFile, datapoint_count, subset_count)
fig.figure.savefig(outFigure)
