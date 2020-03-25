import os, re, csv
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')


'''
Makes a violin plot of the total BDM measures (just BDM, BDM density, or something else), 
grouped by host (or some other group)
'''

bdm_pickle_file = 'bdm_pickles/bdms_proteins_EDSSMat90'

# load in the bdm values
bdms_dict = pickle.load(open(bdm_pickle_file, 'rb'))
grouping = bdm_pickle_file.split('_')[-1]

# reconstruct dict for the plot
for key in bdms_dict.keys():
    del bdms_dict[key]['bdms']

# Complexity for entire string, for all fibers, grouped by host
df = pd.DataFrame(bdms_dict)
df = df.transpose()
df = df.reset_index()
df.columns = ['Tail', 'Group', 'BDM', 'BDM Density']
df['BDM'] = df['BDM'].astype(float)
df['BDM Density'] = df['BDM Density'].astype(float)

fig, ax = plt.subplots()
plt.title('BDM Density Grouped by Host (Grouping = ' + grouping + ')')
fig.set_size_inches(10, 5)
sns.violinplot(x="Group", y="BDM Density", data=df, ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
plt.tight_layout()
plt.savefig('figures_bdm_groups/whole_bdms_density_' + grouping + '.png')
