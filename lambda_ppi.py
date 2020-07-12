from functions import Complexity
import os, re, csv, pickle
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pandas as pd
from scipy import stats
import json


# Get the interactions
# Get the BDM values for each protein in each pair
# get correlation
# plot for fun

# get interactions
interactions_file = os.path.join('virus_host_interactions', 't7')
with open(interactions_file, 'r') as content_file:
    content = content_file.read()
content = re.sub('\n', '', content)
content = re.sub('\'', '"', content)
interactions = json.loads(content)

# Load in the BDM values
virus_bdms = pickle.load(open('bdm_pickles/bdms_proteins_EDSSMat90_data_t7', 'rb'))
host_bdms = pickle.load(open('bdm_pickles/bdms_proteins_EDSSMat90_data_t7_host', 'rb'))


# translate viral protein names to match
viral_bdm_proteins = list(virus_bdms.keys())
for v in viral_bdm_proteins:
    print(v)

v
for key in virus_bdms.keys():
    # this gets 25 / 27 proteins matched, NOT the "spike" one
    new_key = re.sub('(pLVX\-EF1alpha\-|pLXV\-EF1a\-|pLVX\-EF1a\-|nCoV2019\-|IRES\-Puro|2xStrep|nCoV\-2019|\-|\_)', '', key)
    if new_key == '':
        continue
    virus_bdms[new_key] = virus_bdms.pop(key)

# Make the edge tuples
ppi_edges = tuple(zip(viral_proteins, host_proteins))

# calculate all the non-ppis
non_ppi_edges = []
for vnode in viral_proteins:
    for hnode in host_proteins:
        if (vnode, hnode) not in ppi_edges:
            non_ppi_edges.append((vnode, hnode))


# --------- PPI EDGES --------

# for each edge, get the different BDM values
bdm_whole_edges = []
bdm_density_edges = []
bdm_second_edges = []
bdm_third_edges = []

for edge in ppi_edges:

    # virus - host
    # try to match, or else just ignore
    try:
        v_bdm_values = virus_bdms[edge[0]]
        h_bdm_values = host_bdms[edge[1]]
    except:
        continue

    bdm_whole_edges.append((v_bdm_values['whole_bdm'], h_bdm_values['whole_bdm']))
    bdm_density_edges.append((v_bdm_values['bdm_density'], h_bdm_values['bdm_density']))
    bdm_second_edges.append((v_bdm_values['second_order'], h_bdm_values['second_order']))
    bdm_third_edges.append((v_bdm_values['third_order'], h_bdm_values['third_order']))

# turn these into dfs for plotting and stuff
df_whole = pd.DataFrame(list(zip(list(zip(*bdm_whole_edges))[0], list(zip(*bdm_whole_edges))[1], ['Whole']*len(list(zip(*bdm_whole_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])
df_density = pd.DataFrame(list(zip(list(zip(*bdm_density_edges))[0], list(zip(*bdm_density_edges))[1], ['Density']*len(list(zip(*bdm_density_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])
df_second = pd.DataFrame(list(zip(list(zip(*bdm_second_edges))[0], list(zip(*bdm_second_edges))[1], ['Secondary']*len(list(zip(*bdm_second_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])
df_third = pd.DataFrame(list(zip(list(zip(*bdm_third_edges))[0], list(zip(*bdm_third_edges))[1], ['Tertiary']*len(list(zip(*bdm_third_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])

frames = [df_whole, df_density, df_second, df_third]
bdm_df = pd.concat(frames)
bdm_df.reset_index(drop=True, inplace=True)
bdm_df['PPI'] = 1


# --------- NON PPI EDGES --------

# for each edge, get the different BDM values
bdm_whole_edges = []
bdm_density_edges = []
bdm_second_edges = []
bdm_third_edges = []

for edge in non_ppi_edges:

    # virus - host
    # try to match, or else just ignore
    try:
        v_bdm_values = virus_bdms[edge[0]]
        h_bdm_values = host_bdms[edge[1]]
    except:
        continue

    bdm_whole_edges.append((v_bdm_values['whole_bdm'], h_bdm_values['whole_bdm']))
    bdm_density_edges.append((v_bdm_values['bdm_density'], h_bdm_values['bdm_density']))
    bdm_second_edges.append((v_bdm_values['second_order'], h_bdm_values['second_order']))
    bdm_third_edges.append((v_bdm_values['third_order'], h_bdm_values['third_order']))

# turn these into dfs for plotting and stuff
df_whole2 = pd.DataFrame(list(zip(list(zip(*bdm_whole_edges))[0], list(zip(*bdm_whole_edges))[1], ['Whole']*len(list(zip(*bdm_whole_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])
df_density2 = pd.DataFrame(list(zip(list(zip(*bdm_density_edges))[0], list(zip(*bdm_density_edges))[1], ['Density']*len(list(zip(*bdm_density_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])
df_second2 = pd.DataFrame(list(zip(list(zip(*bdm_second_edges))[0], list(zip(*bdm_second_edges))[1], ['Secondary']*len(list(zip(*bdm_second_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])
df_third2 = pd.DataFrame(list(zip(list(zip(*bdm_third_edges))[0], list(zip(*bdm_third_edges))[1], ['Tertiary']*len(list(zip(*bdm_third_edges))[1]))),
                  columns=['Viral Protein BDM', 'Host Protein BDM', 'BDM Type'])

frames = [df_whole2, df_density2, df_second2, df_third2]
bdm_df2 = pd.concat(frames)
bdm_df2.reset_index(drop=True, inplace=True)
bdm_df2['PPI'] = 0


# --------- PUT THEM TOGETHER ---------

# make these into one DF
frames = [bdm_df, bdm_df2]
bdm_df_all = pd.concat(frames)
bdm_df_all = bdm_df_all.sort_values(by='PPI', ascending=True)

# get correlations for all the BDM measures
corr_whole_ppi = stats.spearmanr(df_whole['Viral Protein BDM'], df_whole['Host Protein BDM'])
corr_density_ppi = stats.spearmanr(df_density['Viral Protein BDM'], df_density['Host Protein BDM'])
corr_second_ppi = stats.spearmanr(df_second['Viral Protein BDM'], df_second['Host Protein BDM'])
corr_third_ppi = stats.spearmanr(df_third['Viral Protein BDM'], df_third['Host Protein BDM'])

# and for the non-ppis
corr_whole_non_ppi = stats.spearmanr(df_whole2['Viral Protein BDM'], df_whole2['Host Protein BDM'])
corr_density_non_ppi = stats.spearmanr(df_density2['Viral Protein BDM'], df_density2['Host Protein BDM'])
corr_second_non_ppi = stats.spearmanr(df_second2['Viral Protein BDM'], df_second2['Host Protein BDM'])
corr_third_non_ppi = stats.spearmanr(df_third2['Viral Protein BDM'], df_third2['Host Protein BDM'])

# Do independent t-test for all 4
t_test_whole = stats.ttest_ind(tuple(zip(df_whole['Viral Protein BDM'], df_whole['Host Protein BDM'])),
                tuple(zip(df_whole2['Viral Protein BDM'], df_whole2['Host Protein BDM'])))
t_test_density = stats.ttest_ind(tuple(zip(df_density['Viral Protein BDM'], df_density['Host Protein BDM'])),
                tuple(zip(df_density2['Viral Protein BDM'], df_density2['Host Protein BDM'])))
t_test_second = stats.ttest_ind(tuple(zip(df_second['Viral Protein BDM'], df_second['Host Protein BDM'])),
                tuple(zip(df_second2['Viral Protein BDM'], df_second2['Host Protein BDM'])))
t_test_third = stats.ttest_ind(tuple(zip(df_third['Viral Protein BDM'], df_third['Host Protein BDM'])),
                tuple(zip(df_third2['Viral Protein BDM'], df_third2['Host Protein BDM'])))
t

# plot these x y values and show the correlation value
colors = ["#d6d6d6", "#000000"]
sns.set_palette(sns.color_palette(colors))

g = sns.relplot(x='Viral Protein BDM', y='Host Protein BDM', col='BDM Type', hue='PPI',
                data=bdm_df_all, height=3, kind="scatter", s=1.5, alpha=0.3, edgecolor=None,
                facet_kws={'sharey': False, 'sharex': False}, hue_order=[0, 1])
plt.show()



'''# Make a histogram of the host bdm values
# same with virus

# try to match virus proteins with host ones based on bdm values

# load in the bdm values
bdms_virus = pickle.load(open('bdm_pickles/bdms_proteins_EDSSMat90_data_t7', 'rb'))
bdms_host = pickle.load(open('bdm_pickles/bdms_proteins_EDSSMat90_data_t7_host', 'rb'))

# just get the whole_bdm values
bdms_virus_whole = [i.get('second_order') for i in bdms_virus.values()]
bdms_host_whole = [i.get('second_order') for i in bdms_host.values()]

sns.distplot(bdms_virus_whole)
plt.title('Histogram of Second Order BDM for T7')
plt.xlabel('Second Order BDM')
plt.ylabel('# Occurrences')
#plt.show()

# folder stuff
figure_folder = os.path.join('figures_host_protein')
if not os.path.exists(figure_folder):
    os.makedirs(figure_folder)
plt.savefig(os.path.join(figure_folder, 'hist_t7_second_order.png'))'''
