#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# 23 Jan 2023
# JRA wrote this with the help of ChatGPT. Thanks, ChatGPT bot!


# In[1]:


# Load libraries
import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from matplotlib.patches import Rectangle
import matplotlib.font_manager as font_manager

small_size = 12
medium_size = 16
matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8
plt.rcParams["font.family"] = ["Helvetica Neue", "Times"]
plt.rcParams['font.size'] = medium_size
font_path = '/opt/X11/share/system_fonts/HelveticaNeue.ttc'
font_name = 'Helvetica Neue'
prop = font_manager.FontProperties(fname=font_path)


# In[2]:


# Load data
userInputFile = 'data/mm10_refseq_genes_heatmaps.txt'


inputFile = pd.read_csv(userInputFile, sep='\t')
inputFile

for column_name in inputFile.columns:
    print(column_name)


# In[110]:


# List relevant groups
'''
E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10.bam_RPKM
E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10.bam_RPKM

E14_WT_D2_RNA_D129T02_rep1_trimV5_mm10.bam_RPKM
E14_WT_D2_RNA_D129T10_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D2_RNA_D129T06_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D2_RNA_D129T14_rep2_trimV5_mm10.bam_RPKM

E14_WT_D4_RNA_D129T03_rep1_trimV5_mm10.bam_RPKM
E14_WT_D4_RNA_D129T11_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D4_RNA_D129T07_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D4_RNA_D129T15_rep2_trimV5_mm10.bam_RPKM

E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10.bam_RPKM
E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10.bam_RPKM

E14_D7_DnmtTKO_EzhipWT_RNA_ams78_trimV5_mm10.bam_RPKM
E14_D7_DnmtTKO_EzhipWT_RNA_ams77_trimV5_mm10.bam_RPKM

E14_D7_DnmtTKO_EzhipKO_RNA_ams80_trimV5_mm10.bam_RPKM
E14_D7_DnmtTKO_EzhipKO_RNA_ams79_trimV5_mm10.bam_RPKM
'''


# In[3]:


# Assign relevant data to groups
gene_names = ['name']
D0_WT = ['E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10.bam_RPKM']
D0_TKO = ['E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10.bam_RPKM']
D2_WT = ['E14_WT_D2_RNA_D129T02_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D2_RNA_D129T10_rep2_trimV5_mm10.bam_RPKM']
D2_TKO = ['E14_TKO_D2_RNA_D129T06_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D2_RNA_D129T14_rep2_trimV5_mm10.bam_RPKM']
D4_WT = ['E14_WT_D4_RNA_D129T03_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D4_RNA_D129T11_rep2_trimV5_mm10.bam_RPKM']
D4_TKO = ['E14_TKO_D4_RNA_D129T07_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D4_RNA_D129T15_rep2_trimV5_mm10.bam_RPKM']
D7_WT = ['E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10.bam_RPKM']
D7_TKO = ['E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10.bam_RPKM']
D7_TKO_EzhipWT = ['E14_D7_DnmtTKO_EzhipWT_RNA_ams78_trimV5_mm10.bam_RPKM', 'E14_D7_DnmtTKO_EzhipWT_RNA_ams77_trimV5_mm10.bam_RPKM']
D7_TKO_EzhipKO = ['E14_D7_DnmtTKO_EzhipKO_RNA_ams80_trimV5_mm10.bam_RPKM', 'E14_D7_DnmtTKO_EzhipKO_RNA_ams79_trimV5_mm10.bam_RPKM']


E3_5_ICM = ['BD_E3.5_ICM_RNA_Zhang2018_rep2_GSM2027205.bam_RPKM', 'BD_E3.5_ICM_RNA_Zhang2018_rep1_GSM2027204.bam_RPKM']
E4_0_ICM = ['BD_E4.0_ICM_RNA_Zhang2018_rep2_GSM2027209.bam_RPKM', 'BD_E4.0_ICM_RNA_Zhang2018_rep1_GSM2027208.bam_RPKM']
E5_5_epi = ['BD_E5.5_Epi_RNA_Zhang2018_rep2_GSM2027211.bam_RPKM', 'BD_E5.5_Epi_RNA_Zhang2018_rep1_GSM2027210.bam_RPKM']
E6_5_epi = ['BD_E6.5_Epi_RNA_Zhang2018_rep2_GSM2027215.bam_RPKM', 'BD_E6.5_Epi_RNA_Zhang2018_rep1_GSM2027214.bam_RPKM']

ICM_WT = ['BD_ICM_Dnmt3aWT_RNA_rep1_trimV5_mm10.bam_RPKM', 'BD_ICM_Dnmt3aWT_RNA_rep2_trimV5_mm10.bam_RPKM']
ICM_KO = ['BD_ICM_Dnmt3aKO_RNA_rep1_trimV5_mm10.bam_RPKM', 'BD_ICM_Dnmt3aKO_RNA_rep2_trimV5_mm10.bam_RPKM']

E8_5_epi_1WT = ['C57_E8.5_Dnmt1WT_RNA_Dahlet2020_GSM3752646_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt1WT_RNA_Dahlet2020_GSM3752647_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt1WT_RNA_Dahlet2020_GSM3752648_trimV5_mm10.bam_RPKM']
E8_5_epi_1KO = ['C57_E8.5_Dnmt1KO_RNA_Dahlet2020_GSM3752651_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt1KO_RNA_Dahlet2020_GSM3752652_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt1KO_RNA_Dahlet2020_GSM3752653_trimV5_mm10.bam_RPKM']

Suz12_WT = ['BL6129_2i_mESC_WT_RNA_rep3_Hojfeldt2019_GSM3890926_trimV5_mm10.bam_RPKM', 'BL6129_2i_mESC_WT_RNA_rep2_Hojfeldt2019_GSM3890918_trimV5_mm10.bam_RPKM', 'BL6129_2i_mESC_WT_RNA_rep1_Hojfeldt2019_GSM3890910_trimV5_mm10.bam_RPKM']
Suz12_KO = ['BL6129_2i_mESC_Suz12_KO_RNA_rep3_Hojfeldt2019_GSM3890927_trimV5_mm10.bam_RPKM', 'BL6129_2i_mESC_Suz12_KO_RNA_rep2_Hojfeldt2019_GSM3890919_trimV5_mm10.bam_RPKM', 'BL6129_2i_mESC_Suz12_KO_RNA_rep1_Hojfeldt2019_GSM3890911_trimV5_mm10.bam_RPKM']

Eed_WT = ['E14_mESC_2i_RNA_vanMierlo2019_GSM2711863_trimV5_mm10.bam_RPKM']
Eed_KO = ['E14_mESC_2i_EedKO_RNA_vanMierlo2019_GSM2711865_trimV5_mm10.bam_RPKM']

# replace with ICM Dnmt3aKO!
# E8_5_epi_Het = ['C57_E8.5_Dnmt3aHet_RNA_Dahlet2020_GSM3752654_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3aHet_RNA_Dahlet2020_GSM3752655_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3aHet_RNA_Dahlet2020_GSM3752656_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3aHet_RNA_Dahlet2020_GSM3752657_trimV5_mm10.bam_RPKM']
# E8_5_epi_abdKO = ['C57_E8.5_Dnmt3abdKO_RNA_Dahlet2020_GSM3752658_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3abdKO_RNA_Dahlet2020_GSM3752659_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3abdKO_RNA_Dahlet2020_GSM3752660_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3abdKO_RNA_Dahlet2020_GSM3752661_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3abdKO_RNA_Dahlet2020_GSM3752662_trimV5_mm10.bam_RPKM', 'C57_E8.5_Dnmt3abdKO_RNA_Dahlet2020_GSM3752663_trimV5_mm10.bam_RPKM']


# In[4]:


# Create dataframe with relevant groups
df = inputFile[gene_names + \
               D0_WT + D0_TKO + \
               D2_WT + D2_TKO + \
               D4_WT + D4_TKO + \
               D7_WT + D7_TKO + \
               E3_5_ICM + E4_0_ICM + E5_5_epi + E6_5_epi + \
               ICM_WT + ICM_KO + \
               E8_5_epi_1WT + E8_5_epi_1KO + \
               Suz12_WT + Suz12_KO + \
               Eed_WT + Eed_KO
              ]
df



# In[9]:


### Filter the dataframe to only include the row corresponding to gene X
gene_name = 'Ezh1'
gene = df.loc[df['name'] == gene_name ]

# Get the actual values, not some pandas series thing I don't understand
D0_WT_vals=gene[D0_WT].iloc[0].tolist()
D0_TKO_vals=gene[D0_TKO].iloc[0].tolist()
D2_WT_vals=gene[D2_WT].iloc[0].tolist()
D2_TKO_vals=gene[D2_TKO].iloc[0].tolist()
D4_WT_vals=gene[D4_WT].iloc[0].tolist()
D4_TKO_vals=gene[D4_TKO].iloc[0].tolist()
D7_WT_vals=gene[D7_WT].iloc[0].tolist()
D7_TKO_vals=gene[D7_TKO].iloc[0].tolist()

E3_5_ICM_vals=gene[E3_5_ICM].iloc[0].tolist()
E4_0_ICM_vals=gene[E4_0_ICM].iloc[0].tolist()
E5_5_epi_vals=gene[E5_5_epi].iloc[0].tolist()
E6_5_epi_vals=gene[E6_5_epi].iloc[0].tolist()
E8_5_epi_1WT_vals=gene[E8_5_epi_1WT].iloc[0].tolist()
E8_5_epi_1KO_vals=gene[E8_5_epi_1KO].iloc[0].tolist()
ICM_WT_vals=gene[ICM_WT].iloc[0].tolist()
ICM_KO_vals=gene[ICM_KO].iloc[0].tolist()

Suz12_WT_vals=gene[Suz12_WT].iloc[0].tolist()
Suz12_KO_vals=gene[Suz12_KO].iloc[0].tolist()
Eed_WT_vals=gene[Eed_WT].iloc[0].tolist()
Eed_KO_vals=gene[Eed_KO].iloc[0].tolist()




# Calculate the mean of each group
D0_WT_mean=np.mean(D0_WT_vals)
D0_TKO_mean=np.mean(D0_TKO_vals)
D2_WT_mean=np.mean(D2_WT_vals)
D2_TKO_mean=np.mean(D2_TKO_vals)
D4_WT_mean=np.mean(D4_WT_vals)
D4_TKO_mean=np.mean(D4_TKO_vals)
D7_WT_mean=np.mean(D7_WT_vals)
D7_TKO_mean=np.mean(D7_TKO_vals)

E3_5_ICM_mean=np.mean(E3_5_ICM_vals)
E4_0_ICM_mean=np.mean(E4_0_ICM_vals)
E5_5_epi_mean=np.mean(E5_5_epi_vals)
E6_5_epi_mean=np.mean(E6_5_epi_vals)

E8_5_epi_1WT_mean=np.mean(E8_5_epi_1WT_vals)
E8_5_epi_1KO_mean=np.mean(E8_5_epi_1KO_vals)
ICM_WT_mean=np.mean(ICM_WT_vals)
ICM_KO_mean=np.mean(ICM_KO_vals)

Suz12_WT_mean=np.mean(Suz12_WT_vals)
Suz12_KO_mean=np.mean(Suz12_KO_vals)
Eed_WT_mean=np.mean(Eed_WT_vals)
Eed_KO_mean=np.mean(Eed_KO_vals)




# Define the data and labels
data1 = np.array([[D0_WT_mean, D2_WT_mean, D4_WT_mean, D7_WT_mean], 
                 [D0_TKO_mean, D2_TKO_mean, D4_TKO_mean, D7_TKO_mean]])

data2 = np.array([[ICM_WT_mean, E8_5_epi_1WT_mean],
                  [ICM_KO_mean, E8_5_epi_1KO_mean]] )

data3 = np.array([[E3_5_ICM_mean, E4_0_ICM_mean, E5_5_epi_mean, E6_5_epi_mean]])

data4 = np.array([[Suz12_WT_mean, Eed_WT_mean],
                  [Suz12_KO_mean, Eed_KO_mean]] )

rows = ['WT', 'TKO']
columns1 = ['D0', 'D2', 'D4', 'D7']
columns2 = ['ICM',  'E8.5']
columns3 = ['E3.5', 'E4.0', 'E5.5', ' E6.5']
columns4 = ['Suz12', 'Eed']

# Create the figure and axes objects
fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(15, 5), gridspec_kw={'width_ratios': [len(columns1), len(columns2), len(columns3), len(columns4)]})

# Get max's of each heatmap
vmax1 = np.max(data1) if np.max(data1) > 1 else 1
vmax2 = np.max(data2) if np.max(data2) > 1 else 1
vmax3 = np.max(data3) if np.max(data3) > 1 else 1
vmax4 = np.max(data4) if np.max(data4) > 1 else 1

# Create the heatmaps
im1 = axs[0].imshow(data1, cmap='coolwarm', vmin=0, vmax=vmax1) # bones are their money, and worms are their dollars
im2 = axs[1].imshow(data2, cmap='coolwarm', vmin=0, vmax=vmax2)
im3 = axs[2].imshow(data3, cmap='coolwarm', vmin=0, vmax=vmax3)
im4 = axs[3].imshow(data4, cmap='coolwarm', vmin=0, vmax=vmax4)

# Add the colorbars
cbar1 = axs[0].figure.colorbar(im1, ax=axs[0], shrink=0.295)
cbar2 = axs[1].figure.colorbar(im2, ax=axs[1], shrink=0.295)
cbar3 = axs[2].figure.colorbar(im3, ax=axs[2], shrink=0.15)
cbar4 = axs[3].figure.colorbar(im4, ax=axs[3], shrink=0.295)



# Set the ticks and labels for the first heatmap
axs[0].set_xticks(np.arange(len(columns1)))
axs[0].set_yticks(np.arange(len(rows)))
axs[0].set_xticklabels(columns1, fontname=font_name, fontproperties=prop)
axs[0].set_yticklabels(rows, fontname=font_name, fontproperties=prop)
cbar1.set_ticks([0, round(vmax1, 2)])
plt.setp(cbar1.ax.get_yticklabels(), fontname=font_name, fontproperties=prop, fontsize=small_size)    

# Set the ticks and labels for the second heatmap
axs[1].set_xticks(np.arange(len(columns2)))
axs[1].set_yticks(np.arange(len(rows)))
axs[1].set_xticklabels(columns2, fontname=font_name, fontproperties=prop)
axs[1].set_yticklabels(rows, fontname=font_name, fontproperties=prop)
cbar2.set_ticks([0, round(vmax2, 2)])
plt.setp(cbar2.ax.get_yticklabels(), fontname=font_name, fontproperties=prop, fontsize=small_size)  

# Set the ticks and labels for the third heatmap
axs[2].set_xticks(np.arange(len(columns3)))
axs[2].set_yticks(np.arange(1))
axs[2].set_xticklabels(columns3, fontname=font_name, fontproperties=prop)
cbar3.set_ticks([0, round(vmax3, 2)])
plt.setp(cbar3.ax.get_yticklabels(), fontname=font_name, fontproperties=prop, fontsize=small_size)    

# Set the ticks and labels for the fourth heatmap
axs[3].set_xticks(np.arange(len(columns4)))
axs[3].set_yticks(np.arange(len(rows)))
axs[3].set_xticklabels(columns4, fontname=font_name, fontproperties=prop)
axs[3].set_yticklabels(rows, fontname=font_name, fontproperties=prop)
cbar4.set_ticks([0, round(vmax4, 2)])
plt.setp(cbar1.ax.get_yticklabels(), fontname=font_name, fontproperties=prop, fontsize=small_size) 

# Loop over data dimensions and create text annotations for the first heatmap
max_value_index = np.unravel_index(np.argmax(data1, axis=None), data1.shape)
for i in range(len(data1)):
    for j in range(len(columns1)):
        if (i, j) == max_value_index:
            axs[0].text(j, i, round(data1[i, j], 1), ha='center', va='center', color='black')
        else:
            axs[0].text(j, i, '', ha='center', va='center', color='black')

# Loop over data dimensions and create text annotations for the second heatmap
max_value_index = np.unravel_index(np.argmax(data2, axis=None), data2.shape)
for i in range(len(data2)):
    for j in range(len(columns2)):
        if (i, j) == max_value_index:
            axs[1].text(j, i, round(data2[i, j], 1), ha='center', va='center', color='black')
        else:
            axs[1].text(j, i, '', ha='center', va='center', color='black')

# Loop over data dimensions and create text annotations for the third heatmap
max_value_index = np.unravel_index(np.argmax(data3, axis=None), data3.shape)
for i in range(len(data3)):
    for j in range(len(columns3)):
        if (i, j) == max_value_index:
            axs[2].text(j, i, round(data3[i, j], 1), ha='center', va='center', color='black')
        else:
            axs[2].text(j, i, '', ha='center', va='center', color='black')
            
# Loop over data dimensions and create text annotations for the fourth heatmap
max_value_index = np.unravel_index(np.argmax(data4, axis=None), data4.shape)
for i in range(len(data4)):
    for j in range(len(columns4)):
        if (i, j) == max_value_index:
            axs[3].text(j, i, round(data4[i, j], 1), ha='center', va='center', color='black')
        else:
            axs[3].text(j, i, '', ha='center', va='center', color='black')
            
# Rotate the x-axis labels
plt.setp(axs[0].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(axs[1].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(axs[2].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(axs[3].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
axs[0].set_title(gene_name, fontname=font_name, fontproperties=prop)

outFigure = "/Users/jra/Desktop/%s_heatmap.svg" % (gene_name)
plt.savefig(outFigure)
plt.show()


# Save mean to text file
outText = "/Users/jra/Desktop/%s_meanvals.txt" % (gene_name)
with open(outText, "w") as file:
    file.write("Mean RPKM values"+ "\n")
    file.write("D0 mean: WT: %s \tTKO: %s\n" % (D0_WT_mean, D0_TKO_mean))
    file.write("D2 mean: WT: %s \tTKO: %s\n" % (D2_WT_mean, D2_TKO_mean))
    file.write("D4 mean: WT: %s \tTKO: %s\n" % (D4_WT_mean, D4_TKO_mean))
    file.write("D7 mean: WT: %s \tTKO: %s\n" % (D7_WT_mean, D7_TKO_mean))
#    file.write("D7 mean: TKO_EzhipWT: %s \tTKO_EzhipKO: %s\n" % (D7_TKO_EzhipWT_mean, D7_TKO_EzhipKO_mean))
    file.write("ICM mean: WT: %s \tDnmt3amatKO: %s\n" % (ICM_WT_mean, ICM_KO_mean))
    file.write("E8.5 epi mean: WT: %s \tDnmt1zygKO: %s\n" % (E8_5_epi_1WT_mean, E8_5_epi_1KO_mean))
    file.write("Suz12 mean: WT: %s \tKO: %s\n" % (Suz12_WT_mean, Suz12_KO_mean))
    file.write("Eed mean: WT: %s \tKO: %s\n" % (Eed_WT_mean, Eed_KO_mean))



# In[19]:


# Filter the dataframe to only include the row corresponding to gene X
gene_name = 'Tuba3b'
gene = df.loc[df['name'] == gene_name ]

# Get the actual values, not some pandas series thing I don't understand
D0_WT_vals=gene[D0_WT].iloc[0].tolist()
D0_TKO_vals=gene[D0_TKO].iloc[0].tolist()
D2_WT_vals=gene[D2_WT].iloc[0].tolist()
D2_TKO_vals=gene[D2_TKO].iloc[0].tolist()
D4_WT_vals=gene[D4_WT].iloc[0].tolist()
D4_TKO_vals=gene[D4_TKO].iloc[0].tolist()
D7_WT_vals=gene[D7_WT].iloc[0].tolist()
D7_TKO_vals=gene[D7_TKO].iloc[0].tolist()

E3_5_ICM_vals=gene[E3_5_ICM].iloc[0].tolist()
E4_0_ICM_vals=gene[E4_0_ICM].iloc[0].tolist()
E5_5_epi_vals=gene[E5_5_epi].iloc[0].tolist()
E6_5_epi_vals=gene[E6_5_epi].iloc[0].tolist()
E8_5_epi_1WT_vals=gene[E8_5_epi_1WT].iloc[0].tolist()
E8_5_epi_1KO_vals=gene[E8_5_epi_1KO].iloc[0].tolist()
ICM_WT_vals=gene[ICM_WT].iloc[0].tolist()
ICM_KO_vals=gene[ICM_KO].iloc[0].tolist()

# Calculate the mean of each group
D0_WT_mean=np.mean(D0_WT_vals)
D0_TKO_mean=np.mean(D0_TKO_vals)
D2_WT_mean=np.mean(D2_WT_vals)
D2_TKO_mean=np.mean(D2_TKO_vals)
D4_WT_mean=np.mean(D4_WT_vals)
D4_TKO_mean=np.mean(D4_TKO_vals)
D7_WT_mean=np.mean(D7_WT_vals)
D7_TKO_mean=np.mean(D7_TKO_vals)

E3_5_ICM_mean=np.mean(E3_5_ICM_vals)
E4_0_ICM_mean=np.mean(E4_0_ICM_vals)
E5_5_epi_mean=np.mean(E5_5_epi_vals)
E6_5_epi_mean=np.mean(E6_5_epi_vals)

E8_5_epi_1WT_mean=np.mean(E8_5_epi_1WT_vals)
E8_5_epi_1KO_mean=np.mean(E8_5_epi_1KO_vals)
ICM_WT_mean=np.mean(ICM_WT_vals)
ICM_KO_mean=np.mean(ICM_KO_vals)

# Define the data and labels
data1 = np.array([[D0_WT_mean, D2_WT_mean, D4_WT_mean, D7_WT_mean]])

data2 = np.array([[ICM_WT_mean, E8_5_epi_1WT_mean]] )

data3 = np.array([[E3_5_ICM_mean, E4_0_ICM_mean, E5_5_epi_mean, E6_5_epi_mean]])

rows = ['WT']
columns1 = ['D0', 'D2', 'D4', 'D7']
columns2 = ['ICM',  'E8.5']
columns3 = ['E3.5', 'E4.0', 'E5.5', ' E6.5']

# Create the figure and axes objects
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 5), gridspec_kw={'width_ratios': [len(columns1), len(columns2), len(columns3)]})

# Get max's of each heatmap
vmax1 = np.max(data1) if np.max(data1) > 1 else 1
vmax2 = np.max(data2) if np.max(data2) > 1 else 1
vmax3 = np.max(data3) if np.max(data3) > 1 else 1


# Create the heatmaps
im1 = axs[0].imshow(data1, cmap='coolwarm', vmin=0, vmax=vmax1) # bones are their money, and worms are their dollars
im2 = axs[1].imshow(data2, cmap='coolwarm', vmin=0, vmax=vmax2)
im3 = axs[2].imshow(data3, cmap='coolwarm', vmin=0, vmax=vmax3)

# Add the colorbars
cbar1 = axs[0].figure.colorbar(im1, ax=axs[0], shrink=0.295)
cbar2 = axs[1].figure.colorbar(im2, ax=axs[1], shrink=0.295)
cbar3 = axs[2].figure.colorbar(im3, ax=axs[2], shrink=0.15)




# Set the ticks and labels for the first heatmap
axs[0].set_xticks(np.arange(len(columns1)))
axs[0].set_yticks(np.arange(len(rows)))
axs[0].set_xticklabels(columns1, fontname=font_name, fontproperties=prop)
axs[0].set_yticklabels(rows, fontname=font_name, fontproperties=prop)
cbar1.set_ticks([0, round(vmax1, 2)])
plt.setp(cbar1.ax.get_yticklabels(), fontname=font_name, fontproperties=prop, fontsize=small_size)    

# Set the ticks and labels for the second heatmap
axs[1].set_xticks(np.arange(len(columns2)))
axs[1].set_yticks(np.arange(len(rows)))
axs[1].set_xticklabels(columns2, fontname=font_name, fontproperties=prop)
axs[1].set_yticklabels(rows, fontname=font_name, fontproperties=prop)
cbar2.set_ticks([0, round(vmax2, 2)])
plt.setp(cbar2.ax.get_yticklabels(), fontname=font_name, fontproperties=prop, fontsize=small_size)  

# Set the ticks and labels for the third heatmap
axs[2].set_xticks(np.arange(len(columns3)))
axs[2].set_yticks(np.arange(1))
axs[2].set_xticklabels(columns3, fontname=font_name, fontproperties=prop)
cbar3.set_ticks([0, round(vmax3, 2)])
plt.setp(cbar3.ax.get_yticklabels(), fontname=font_name, fontproperties=prop, fontsize=small_size)    

# Loop over data dimensions and create text annotations for the first heatmap
max_value_index = np.unravel_index(np.argmax(data1, axis=None), data1.shape)
for i in range(len(data1)):
    for j in range(len(columns1)):
            axs[0].text(j, i, round(data1[i, j], 1), ha='center', va='center', color='black')


# Loop over data dimensions and create text annotations for the second heatmap
max_value_index = np.unravel_index(np.argmax(data2, axis=None), data2.shape)
for i in range(len(data2)):
    for j in range(len(columns2)):
            axs[1].text(j, i, round(data2[i, j], 1), ha='center', va='center', color='black')

# Loop over data dimensions and create text annotations for the third heatmap
max_value_index = np.unravel_index(np.argmax(data3, axis=None), data3.shape)
for i in range(len(data3)):
    for j in range(len(columns3)):
            axs[2].text(j, i, round(data3[i, j], 1), ha='center', va='center', color='black')

            
# Rotate the x-axis labels
plt.setp(axs[0].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(axs[1].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(axs[2].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
axs[0].set_title(gene_name, fontname=font_name, fontproperties=prop)

outFigure = "/Users/jra/Desktop/%s_heatmap_noKO.svg" % (gene_name)
plt.savefig(outFigure)
plt.show()


