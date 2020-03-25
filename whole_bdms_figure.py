import os, re, csv
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')


files = os.listdir('.')
files = list(filter(lambda x: re.search('bdms_proteins', x), files))
files = ['pickle_jar/bdms_proteins_EDSSMat90']

for file in files:

    # load in the bdm values
    bdms_dict = pickle.load(open(file, 'rb'))
    grouping = file.split('_')[-1]

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
    plt.savefig('figures2/whole_bdms_density_' + grouping + '.png')
