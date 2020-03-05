import os, re, csv
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')


files = os.listdir('.')
files = list(filter(lambda x: re.search('bdms_proteins', x), files))

for file in files:

    # load in the bdm values
    bdms_dict = pickle.load(open(file, 'rb'))
    grouping = file.split('_')[2]
    window_size = file.split('_')[3]

    # reconstruct dict for the plot
    for key in bdms_dict.keys():
        del bdms_dict[key]['bdms']

    # Complexity for entire string, for all fibers, grouped by host
    df = pd.DataFrame(bdms_dict)
    df = df.transpose()
    df = df.reset_index()
    df.columns = ['Tail', 'Group', 'BDM']
    df['BDM'] = df['BDM'].astype(float)

    fig, ax = plt.subplots()
    plt.title('BDM Grouped by Host (' + grouping + ', window = ' + window_size + ')')
    fig.set_size_inches(10, 5)
    sns.violinplot(x="Group", y="BDM", data=df, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
    plt.tight_layout()
    plt.savefig('figures2/whole_bdms_' + grouping + '_' + window_size + '.png')
