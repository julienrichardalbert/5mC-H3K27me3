#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None # stop annoying warning when clipping columns. default='warn'

# normally, I install python packages with conda while in the "jupyter" conda environment I made
# however, combat isn't packaged in conda, so I do this:
#### !{sys.executable} -m pip install combat
# as recommended by:
# https://jakevdp.github.io/blog/2017/12/05/installing-python-packages-from-jupyter/
# this will install combat. I thought I had to do this each time I launch a notebook but I do not! Yay!!
# but combat (pypi) is listed in my conda environment.... I'm confused. Do pip installs show up in conda list?
from combat.pycombat import pycombat


# In[2]:


# prepare data
# the indexes correspond to the gene names
# the column names correspond to the sample names
df_expression = pd.read_csv("data/mm10_refGene_PCA.txt", sep='\t', index_col=0)

print(df_expression.shape)
#print(df_expression.columns)
#print(df_expression.head(2))


# In[3]:


# assign batches manually
batch = [ 1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4]
# make colours
colour_names = {0:'red', 1:'blue', 2:'grey', 3:'yellow', 4:'green', 5:'purple', 6:'pink', 7:'orange', 8:'teal' }

colour_indexes = []
for x in batch:
    colour_indexes.append(colour_names.get(x))
print(colour_indexes)

# get 'labels' for the x-axis
labels = list(df_expression.columns)
print(labels)
print(len(batch))


# In[4]:


df_expression = df_expression.loc[(df_expression!=0).any(1)] # drop rows that have all 0 values
print("Dropping rows with only 0 values")
print(df_expression.shape)


#df_expression
# SOMETHING IS WRONG WITH THE TABLE LETS FIND OUT WHAT
df_expression.isnull().values.any()
df_expression.isnull().sum()
null_data = df_expression[df_expression.isnull().any(axis=1)]
null_data
df_expression


# In[5]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (20,4)


boxPlot = plt.boxplot(df_expression, \
                    labels=labels, \
                    patch_artist=True, notch=True, \
                    medianprops={'color':'Black'}, \
                    flierprops={'markersize':2, 'marker':'o', 'markerfacecolor':'white', 'markeredgecolor':'black', 'markeredgewidth':0.5})
plt.xticks(rotation=90)
for patch, color in zip(boxPlot['boxes'], colour_indexes):
    patch.set_facecolor(color)
plt.ylim(-2.1, None)    

 
#outFigure = "/Users/jra/Desktop/bigpapa_raw.svg"
#plt.savefig(outFigure)
plt.show() # showing the figure "moves" the figure to Shell, leaving nothing behind. Show after saving file.


# In[6]:


df_log = np.log2(df_expression.clip(upper=None, lower=0)+1)
boxPlot = plt.boxplot(df_log, \
                    labels=labels, \
                    patch_artist=True, notch=True, \
                    medianprops={'color':'Black'}, \
                    flierprops={'markersize':2, 'marker':'o', 'markerfacecolor':'white', 'markeredgecolor':'black', 'markeredgewidth':0.5})
plt.xticks(rotation=90)
for patch, color in zip(boxPlot['boxes'], colour_indexes):
    patch.set_facecolor(color)
plt.ylim(-2.1, None)    

 
#outFigure = "/Users/jra/Desktop/bigpapa_raw.svg"
#plt.savefig(outFigure)
plt.show() # showing the figure "moves" the figure to Shell, leaving nothing behind. Show after saving file.


# In[7]:


# run pyComBat
df_corrected = pycombat(df_expression, batch)
#df_corrected = pycombat(df_expression, batch)

df_corrected


# In[8]:


boxPlot = plt.boxplot(df_corrected, \
                    labels=labels, \
                    patch_artist=True, notch=True, \
                    medianprops={'color':'Black'}, \
                    flierprops={'markersize':2, 'marker':'o', 'markerfacecolor':'white', 'markeredgecolor':'black', 'markeredgewidth':0.5})
plt.xticks(rotation=90)
for patch, color in zip(boxPlot['boxes'], colour_indexes):
    patch.set_facecolor(color)
plt.ylim(-2.1, None)    
 
outFigure = "/Users/jra/Desktop/bigpapa_raw.svg"
plt.savefig(outFigure)
plt.show() # showing the figure "moves" the figure to Shell, leaving nothing behind. Show after saving file.


# In[9]:


df_corrected_log = np.log2(df_corrected.clip(upper=None, lower=0)+1)
boxPlot = plt.boxplot(df_corrected_log, \
                    labels=labels, \
                    patch_artist=True, notch=True, \
                    medianprops={'color':'Black'}, \
                    flierprops={'markersize':2, 'marker':'o', 'markerfacecolor':'white', 'markeredgecolor':'black', 'markeredgewidth':0.5})
plt.xticks(rotation=90)
for patch, color in zip(boxPlot['boxes'], colour_indexes):
    patch.set_facecolor(color)
plt.ylim(-2.1, None)    
 
#outFigure = "/Users/jra/Desktop/bigpapa_raw.svg"
#plt.savefig(outFigure)
plt.show() # showing the figure "moves" the figure to Shell, leaving nothing behind. Show after saving file.


# In[12]:


np.savez_compressed('/Users/jra/Desktop/vivo_vitro_RNA_combatCorrected.npz', matrix=df_corrected, labels=labels)
np.savez_compressed('/Users/jra/Desktop/vivo_vitro_RNA_scalp_raw.npz', matrix=df_expression, labels=labels)
pd.DataFrame(df_corrected).to_csv("/Users/jra/Desktop/vivo_vitro_RNA_combatCorrected.txt", sep="\t")
pd.DataFrame(df_expression).to_csv("/Users/jra/Desktop/vivo_vitro_RNA_scalp_raw.txt", sep="\t")


# In[10]:


from sklearn.decomposition import PCA
plt.rcParams["figure.figsize"] = (12,12)

samples = df_expression.transpose()
pca = PCA(2)  # project from 64 to 2 dimensions
projected = pca.fit_transform(samples)
print(df_expression.shape)
print(samples.shape)
print(projected.shape)



plt.scatter(projected[:, 0], projected[:, 1],
            c=batch, edgecolor='none', alpha=1,
            cmap=plt.cm.get_cmap('tab20b', 8), s=200)

x_title = 'PC1: ' + (pca.explained_variance_ratio_ * 100)[0].round(2).astype(str) + '% of variation explained'
y_title = 'PC2: ' + (pca.explained_variance_ratio_ * 100)[1].round(2).astype(str) + '% of variation explained'
plt.xlabel(x_title)
plt.ylabel(y_title)
plt.title('made by Julien using python in a jupyter notebook')

for i, label in enumerate(labels):
    plt.annotate(label, (projected[:,0][i], projected[:,1][i]))
outFigure = "/Users/jra/Desktop/pythonPCA_raw2.svg"
plt.savefig(outFigure)


# In[11]:


plt.rcParams["figure.figsize"] = (12,12)

samples = df_corrected.transpose()
pca = PCA(2)  # project from 64 to 2 dimensions
projected = pca.fit_transform(samples)
print(df_corrected.shape)
print(samples.shape)
print(projected.shape)


plt.scatter(projected[:, 0], projected[:, 1],
            c=batch, edgecolor='none', alpha=1,
            cmap=plt.cm.get_cmap('tab20b', 8), s=200)

x_title = 'PC1: ' + (pca.explained_variance_ratio_ * 100)[0].round(2).astype(str) + '% of variation explained'
y_title = 'PC2: ' + (pca.explained_variance_ratio_ * 100)[1].round(2).astype(str) + '% of variation explained'
plt.xlabel(x_title)
plt.ylabel(y_title)
plt.title('made by Julien using python in a jupyter notebook')
for i, label in enumerate(labels):
    plt.annotate(label, (projected[:,0][i], projected[:,1][i]))
#plt.colorbar()
outFigure = "/Users/jra/Desktop/pythonPCA_combat2.svg"
plt.savefig(outFigure)


# In[12]:


# https://stackoverflow.com/questions/22984335/recovering-features-names-of-explained-variance-ratio-in-pca-with-sklearn
# Redo the PCA to find the genes that contribute most variability to our the two PCs
samples = df_corrected.transpose()
pca = PCA(2)  # project from 64 to 2 dimensions
projected = pca.fit_transform(samples)
print(df_corrected.shape)
print(samples.shape)
print(projected.shape)


componentinos = pd.DataFrame(pca.components_, index = ['PC1','PC2'], columns = df_corrected.index)
print(componentinos)
comp_abs = componentinos.transpose().abs()
print(comp_abs)
print(comp_abs[['PC1']].idxmax())
print(comp_abs[['PC2']].idxmax())
comp_abs.sort_values(by=['PC1'], ascending=False)


# In[56]:


np.savez_compressed('/Users/jra/Desktop/vivo_vitro_RNA__PCcomponents.npz', matrix=comp_abs, labels=labels)
pd.DataFrame(comp_abs).to_csv("/Users/jra/Desktop/vivo_vitro_RNA_PCcomponents.txt", sep="\t")

