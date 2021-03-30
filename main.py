# TODO: Need to make sure all the directories self-create
# TODO: Need to make sure all dependencies can auto-install, cut out any unnecessary ones
# TODO: Make the main script, so it can run from the command line
# TODO: set up the first-time calculation stuff

# assume data is downloaded in data folder
# assume PPIs are in data folder too
# For each protein, calculate the BDM (just do whole BDM for now)

import os
import pickle, re
import pandas as pd
import dask.dataframe as dd
import random, json, csv, time
from ete3 import NCBITaxa
import networkx as nx
import numpy as np
from matplotlib.pyplot import figure
from multiprocessing import Pool, freeze_support
import requests

from functions import Complexity

import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('ggplot')
#sns.set(rc={'figure.figsize': (8, 6)})
from networkx.drawing.nx_agraph import graphviz_layout



# =================== Load in all the protein names from the bdm df ===================  # TODO: Need to change the order of loading in files

data_source = 'virusstring'
pickle_out = os.path.join('pickle_jar/', data_source + '_bdms')
"""with open(os.path.join(pickle_out, 'df.p'), 'rb') as f:
    df = pickle.load(f)

# just get the virus protein names
df = df.loc[df['protein_type'] == 'virus']
virus_proteins = list(df['string_id'])


# =================== Get the PPIs from the API ===================

api_call = "https://version-11-0.string-db.org/api/json/interaction_partners?identifiers="

ppis = []

# do ten per call to save time
def chunks(lst, n):
    "Yield successive n-sized chunks from lst."
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

virus_proteins = chunks(virus_proteins, 20)

for proteins in virus_proteins:
    api_call_end = '%0d'.join(proteins)
    r = requests.get(api_call+api_call_end)
    interactions = json.loads(r.text)
    ppis.extend(interactions)
    time.sleep(1)
    print('.')

# save to a pickle file
with open(os.path.join(pickle_out, 'ppis.p'), 'wb') as handle:
    pickle.dump(ppis, handle)

quit()"""

with open(os.path.join(pickle_out, 'ppis.p'), 'rb') as f:
    ppis = pickle.load(f)

# Look at the complexity of these proteins for binning purposes real quick
# initialize class from Functions.py file (only used to calculate bdms)
type = 'proteins'
grouping = 'EDSSMat90'
protein_sequences = 'data'

"""
#genecomplexity = Complexity()
#bdm = genecomplexity.init_bdm(type=type, grouping=grouping)
pickle_out = 'bins_out'

# run the calculation (takes hours, saves result to pickle files, only need to one once per set of proteins)
#genecomplexity.calculate_whole_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory=protein_sequences)
#quit()

# turn into df for plots
with open('bins_out/whole_bdm.p', 'rb') as f:
    bdms = pickle.load(f)

for key in bdms.keys():
    bin = '_'.join(key.split('_')[:-1])
    bdms[key]['bin'] = bin

df = pd.DataFrame.from_dict(bdms, orient='index')
df = df.reset_index(drop=True)

# make the plot
sns.set(font_scale=0.75)
ax = sns.violinplot(x="bin", y="whole_bdm", data=df, scale="width", linewidth=0.2)
[tick.set_rotation(-90) for tick in ax.get_xticklabels()]
plt.tight_layout()
plt.show()
"""



# run metabolic
#perl /slowdata/data1/Genome_profile_software3/METABOLIC-G.pl -in-gn /slowdata/data1/Genome_profile_software3/Rifle_genomes/drep_outout_85/dereplicated_genomes -t 40 -o test_out -m /slowdata/data1/Genome_profile_software

# specify file paths

type = 'proteins'
grouping = 'EDSSMat90'
data_source = 'virusstring'
#ppis = os.path.join('data', data_source + '_ppis')
protein_sequences = os.path.join('data/', data_source + '_sequences')
pickle_out = os.path.join('pickle_jar/', data_source + '_bdms')


# ==================================================================================================
# --------------------------- Step 0: Load in ids, ppis, and taxonomy tables ---------------------------
# ==================================================================================================

# load in the taxonomy tree from ncbi
# When matching, don't go lower than the family level

print('step 0: Load in files')

# load in ID translation tables from string to ncbi
#id_table = dd.read_csv('all_organisms.name_2_string.tsv', delimiter='\t', skiprows=1, header=None, dtype='str')
#id_table = id_table.set_index(2)
#id_table = id_table.compute()
#id_table = id_table.to_dict()

# NCBI taxonomy and tree information
#ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()  # Only for the first time running

# make PPI out-edge dict
# TODO: refactor into other file
def make_ppi_out_edges():

    ppi_out_edges = {}
    with open(os.path.join(ppis, 'protein.links.full.v10.5.txt'), 'rb') as f:

        # skip first line, save header
        header = next(f)
        header = header.split()
        n = 0

        # read in line by line
        for i, line in enumerate(f):

            line = line.split()

            # only load in the experiment and database ones
            experiments = int(line[9])
            database = int(line[11])
            if experiments == 0 and database == 0:
                continue

            p1 = line[0].decode("utf-8")
            p2 = line[1].decode("utf-8")

            # if new key, add as list
            try:
                ppi_out_edges[p1].append(p2)
                ppi_out_edges[p1] = list(set(ppi_out_edges[p1]))
            except:
                ppi_out_edges[p1] = [p2]

            #if p1_ncbi not in ppi_out_edges.keys():
            #    ppi_out_edges[p1_ncbi] = [p2_ncbi]
            #else:  # else append to existing list
            #    ppi_out_edges[p1_ncbi].append(p2_ncbi)

            # bi-directional edges, so add to both keys
            try:
                ppi_out_edges[p2].append(p1)
                ppi_out_edges[p2] = list(set(ppi_out_edges[p2]))
            except:
                ppi_out_edges[p2] = [p1]

            #if p2_ncbi not in ppi_out_edges.keys():
            #    ppi_out_edges[p2_ncbi] = [p1_ncbi]
            #else:  # else append to existing list
            #    ppi_out_edges[p2_ncbi].append(p1_ncbi)

    # pickle this dict
    folder = os.path.join(pickle_out)
    with open(os.path.join(pickle_out, 'ppi_out_edges.p'), 'wb') as handle:
        pickle.dump(ppi_out_edges, handle)


def make_ppi_out_edges2():

    ppi_out_edges = {}
    for ppi in p:

        # skip first line, save header
        header = next(f)
        header = header.split()
        n = 0

        # read in line by line
        for i, line in enumerate(f):

            line = line.split()

            # only load in the experiment and database ones
            experiments = int(line[9])
            database = int(line[11])
            if experiments == 0 and database == 0:
                continue

            p1 = line[0].decode("utf-8")
            p2 = line[1].decode("utf-8")

            # if new key, add as list
            try:
                ppi_out_edges[p1].append(p2)
                ppi_out_edges[p1] = list(set(ppi_out_edges[p1]))
            except:
                ppi_out_edges[p1] = [p2]

            #if p1_ncbi not in ppi_out_edges.keys():
            #    ppi_out_edges[p1_ncbi] = [p2_ncbi]
            #else:  # else append to existing list
            #    ppi_out_edges[p1_ncbi].append(p2_ncbi)

            # bi-directional edges, so add to both keys
            try:
                ppi_out_edges[p2].append(p1)
                ppi_out_edges[p2] = list(set(ppi_out_edges[p2]))
            except:
                ppi_out_edges[p2] = [p1]

            #if p2_ncbi not in ppi_out_edges.keys():
            #    ppi_out_edges[p2_ncbi] = [p1_ncbi]
            #else:  # else append to existing list
            #    ppi_out_edges[p2_ncbi].append(p1_ncbi)

    # pickle this dict
    folder = os.path.join(pickle_out)
    with open(os.path.join(pickle_out, 'ppi_out_edges.p'), 'wb') as handle:
        pickle.dump(ppi_out_edges, handle)


#make_ppi_out_edges2()
#quit()

#with open(os.path.join(pickle_out, 'ppi_out_edges.p'), 'rb') as f:
#    ppi_out_edges = pickle.load(f)

# ==================================================================================================
# --------------------------- Step 1: Get the BDMs of the proteins ---------------------------
# ==================================================================================================

print('step 1: Calculate or load in BDM values')

# initialize class from Functions.py file (only used to calculate bdms)
#genecomplexity = Complexity()
#bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

# run the calculation (takes hours, saves result to pickle files, only need to one once per set of proteins)
#genecomplexity.calculate_whole_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory=protein_sequences)
#quit()

# or load in bdm values from pickle files from previous run
files = os.listdir(pickle_out)
files = list(filter(lambda x: not re.search('DS_Store', x) and not re.search('diffs', x)
                               and re.search('whole_bdm', x), files)) # <<<<<<<<<<<<

# ==================================================================================================
# --------------------------- Step 2: Make a BDMs dict/DB ---------------------------
# ==================================================================================================

print('step 2: Make df of values for each protein')

# Save BDM values to a dict, which gets turned into a DF for easy plotting in seaborn
# load in the results from metabolic amd add these to the df

def make_df():

    dict_to_df = {}

    for file in files:

        print(file)  # <<<<<<<<<<<<<<<<<dev>>>>>>>>>>>>>

        with open(os.path.join(pickle_out, file), 'rb') as f:

            bdms = pickle.load(f)

            for k in bdms.keys():

                # get NCBI id
                if k in id_table[0].keys():
                    ncbi_id = id_table[0][k]  # the first one is the NCBI number
                else:
                    ncbi_id = None

                protein_type = None
                species = None
                family = None
                order = None
                clas = None
                phylum = None
                kingdom = None
                out_edges = None

                # only interested in ones with known names, but keep non-matches in dict anyways
                if ncbi_id:

                    # get species name
                    species = ncbi.get_taxid_translator([int(ncbi_id)])
                    species = list(species.values())[0]

                    # decide if virus or host
                    if species:
                        if re.search('virus', species) or re.search('phag', species) or re.search('vira', species):
                            protein_type = 'virus'
                        else:
                            protein_type = 'host'
                    else:
                        if re.search('virus', bdms[k]['group']) or re.search('phag', bdms[k]['group']) or re.search(
                                'vira', bdms[k]['group']):
                            protein_type = 'virus'
                        else:
                            protein_type = 'host'

                    # get taxonomy information of host
                    if protein_type == 'host':
                        desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                        ranks = Complexity.get_desired_ranks(ncbi_id, desired_ranks)
                        family = ranks['family_id']
                        order = ranks['order_id']
                        clas = ranks['class_id']
                        phylum = ranks['phylum_id']
                        kingdom = ranks['kingdom_id']

                    # if virus, search virus taxonomy file
                    #if protein_type == 'virus':

                        # try to match to a row
                        #for i, r in virus_taxonomy_tree.iterrows():

                            #bdms[k]['family'] = virus_taxonomy_tree[virus_taxonomy_tree['Species'] == species_id].index
                            #bdms[k]['order'] = bdms
                            #bdms[k]['class'] = bdms
                            #bdms[k]['phylum'] = bdms
                            #bdms[k]['kingdom'] = bdms

                else:

                    # decide if virus or host
                    if re.search('virus', bdms[k]['group']) or re.search('phag', bdms[k]['group']) or re.search('vira', bdms[k]['group']):
                        protein_type = 'virus'
                    else:
                        protein_type = 'host'

                    ncbi_id = None
                    species = None
                    family = None
                    order = None
                    clas = None
                    phylum = None
                    kingdom = None

                # get out-edges of PPI network (just ncbi id)
                try:
                    out_edges = ppi_out_edges[k]
                except:
                    out_edges = None

                # save all this information to dict_to_df
                dict_to_df[k] = {}
                dict_to_df[k]['string_id'] = k
                dict_to_df[k]['group'] = bdms[k]['group']
                dict_to_df[k]['whole_bdm'] = bdms[k]['whole_bdm']
                dict_to_df[k]['length'] = bdms[k]['length']
                dict_to_df[k]['protein_type'] = protein_type
                dict_to_df[k]['ncbi_id'] = ncbi_id
                dict_to_df[k]['species'] = species
                dict_to_df[k]['family'] = family
                dict_to_df[k]['order'] = order
                dict_to_df[k]['class'] = clas
                dict_to_df[k]['phylum'] = phylum
                dict_to_df[k]['kingdom'] = kingdom
                dict_to_df[k]['out_edges'] = out_edges

    df = pd.DataFrame.from_dict(dict_to_df, orient='index')
    df = df.reset_index(drop=True)

    # read in metabolic output
    metabolic = pd.read_csv('data/output/METABOLIC_result_each_spreadsheet/METABOLIC_result_worksheet1.tsv', sep='\t')
    df_to_join = {}

    # make new columns
    #df["category"] = ""
    #df["function"] = ""

    # add in the values
    # for each cell in the last column, loop over the values and add category and function to graph
    for index, row in metabolic.iterrows():
        category = row[0]
        function = row[1]
        protein_ids = row[-1].split(',')

        for protein in protein_ids:
            df_to_join[protein] = {}
            df_to_join[protein]['string_id'] = protein
            df_to_join[protein]['category'] = category
            df_to_join[protein]['function'] = function

    df_to_join = pd.DataFrame.from_dict(df_to_join, orient='index')

    # do a join on the two
    df = df.set_index('string_id').join(df_to_join.set_index('string_id'))
    df = df.reset_index()

    # pickle this DF to load in all at once (must have lots of ram)
    with open(os.path.join(pickle_out, 'df.p'), 'wb') as handle:
        pickle.dump(df, handle)

make_df()
#quit()
# TODO: Make all these pickle loading/unloading into a log file that sees if these files are all made or not
#  and whether or not everything needs to be calculated or not

#with open(os.path.join(pickle_out, 'df.p'), 'rb') as f:
#    df = pickle.load(f)


# ==================================================================================================
# --------------------------- Step 3: Make entire PPI network with node attributes -----------------
# ==================================================================================================

print('step 3: Make the entire PPI network')

def make_graph(df):

    # -----> Make the PPI networks

    # make new df for the graph, has to be one column per out node
    node_id = []
    out_node_id = []

    for index, row in df.iterrows():
        if row['out_edges'] == None:
            continue
        for out_node in row['out_edges']:
            node_id.append(row['string_id'])
            out_node_id.append(out_node)

    for_df_graph = list(zip(node_id, out_node_id))
    df_graph = pd.DataFrame(for_df_graph, columns=['node_id', 'out_node_id'])
    graph = nx.from_pandas_edgelist(df_graph, source='node_id', target='out_node_id')

    # add in node attributes from the df
    attributes = df.set_index('string_id').to_dict('index')
    nx.set_node_attributes(graph, attributes)

    # BDM vs node degree, color nodes by protein type
    degrees = []
    graph_degree_dict = dict(graph.degree())
    for protein in df['string_id']:
        try:
            if graph_degree_dict[protein]:
                degrees.append(graph.degree(protein))
        except:
            degrees.append(0)

    df['degree'] = degrees
    # can add more network features to this

    # pickle this graph
    with open(os.path.join(pickle_out, 'graph.p'), 'wb') as handle:
        pickle.dump(graph, handle)

    # update that df
    with open(os.path.join(pickle_out, 'df.p'), 'wb') as handle:
        pickle.dump(df, handle)

#make_graph(df)
#quit()

#with open(os.path.join(pickle_out, 'graph.p'), 'rb') as f:
#    graph = pickle.load(f)

with open(os.path.join(pickle_out, 'df.p'), 'rb') as f:
    df = pickle.load(f)

# 11320.NRAM_I34A1 virus protein for sanity checks

# now i need to check which nodes are missing from the df
#missing_proteins = []
#df_nodes = dict(zip(list(df['string_id']), list(df['string_id'])))
#for node in list(graph.nodes()):
#    try:
#        if df_nodes[node]:
#            continue
#    except:
#        missing_proteins.append(node)

#with open('missing_proteins.txt', 'w') as f:
#    for node in missing_proteins:
#        f.write("%s\n" % node)
#quit()


# ==================================================================================================
# --------------------------- Step 4: Make plots ---------------------------
# ==================================================================================================

print('step 4: Making plots')
# logscale complexity ok!

# TODO: Viruses in the DB aren't in the PPI edge file!!! Even though they appear to have interactions on the website

# to help with dev
#sample = df.sample(10000)
plt.clf()


# ========= Histograms of BDMs =================

# # -----> both the hosts and the viruses
#ax = sns.histplot(data=df, x="whole_bdm", hue="protein_type", stat="density", common_norm=False, element="step", log_scale=(True, False))
#ax.set(xlabel='C', ylabel='Frequency')
#plt.legend(title="Protein Type")
#plt.tight_layout()
#plt.show()
#plt.savefig('bdm_hists.pdf')
#plt.clf()

# -----> both the hosts and the viruses, but stacked
# ax = sns.kdeplot(data=df, x="whole_bdm", hue="protein_type", multiple="fill", log_scale=(True, True))
# ax.set(xlabel='BDM Value', ylabel='Density Fraction', title='Viral + Host Proteins')
# plt.savefig('Figures/both_bdm_hists_fraction.pdf')
# plt.clf()

# # -----> just the viruses, color by group
#sns.set(rc={'figure.figsize': (10, 7)})
#df = df.loc[df['protein_type'] == 'virus']
#ax = sns.displot(data=df, x="whole_bdm", hue="group",
#    log_scale=(True, False), kind="kde", multiple="fill",
#    legend=False, linewidth=0.1, palette="pastel", aspect=2)
#ax.set(xlabel='C', ylabel='Frequency')
#plt.tight_layout()
#plt.show()
#plt.savefig('virus_bdm_hist_group.pdf')
#plt.clf()

# -----> just the hosts, but color by taxonomy
#sns.set(rc={'figure.figsize': (10, 7)})
#df = df.loc[df['protein_type'] == 'host']
#df = df.dropna(subset=['class'])  # get rid of nans for interested column
#df['class'] = df['class'].astype(int)  # make tax_ids integers TODO: for whole df
#sns.displot(data=df, x="whole_bdm", hue="class",
#    log_scale=(True, False), kind="kde", multiple="fill", legend=False,
#    linewidth=0.1, palette="pastel", aspect=1.5)
#plt.xlabel('C')
#plt.ylabel('Frequency')
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=2)
#plt.tight_layout()
#plt.show()
#plt.savefig('host_bdm_hist_class.pdf')
#plt.clf()


# ========= BDM vs Metabolic function ===============

# function (or category) vs bdm value violin plots, make this one cooler!
#df = df.dropna()
#df = df.loc[df.whole_bdm < 10000]  # remove outliers
#sns.set(font_scale=0.5)
#ax = sns.violinplot(x="category", y="whole_bdm", data=df, scale="width", linewidth=0.2)
#[tick.set_rotation(-90) for tick in ax.get_xticklabels()]
#plt.tight_layout()
#plt.xlabel("Metabolic Category")
#plt.ylabel("C")
#plt.show()
#plt.savefig("category_bdm.pdf")

# just for viruses
#sns.set(rc={'figure.figsize': (6, 10)})
#df = df.loc[df['protein_type'] == 'virus']  # 74K proteins here,
#df = df.dropna(subset=['function'])  # but only 889 are annotated
#sns.set(font_scale=0.55)
#ax = sns.violinplot(x="function", y="whole_bdm", data=df, scale="width", linewidth=0.2)
#[tick.set_rotation(-90) for tick in ax.get_xticklabels()]
#plt.tight_layout()
#plt.xlabel("Metabolic Function")
#plt.ylabel("C")
#plt.show()
#plt.savefig("function_bdm_viruses.pdf")
#quit()

# show hosts and viruses at the same time
#df = df.loc[df.whole_bdm < 20000]  # remove outliers
#sns.set(rc={'figure.figsize': (5, 15)})
#sns.set(font_scale=0.5)
#plt.figure(figsize=(8, 12))
#ax = sns.violinplot(x="whole_bdm", y="function", hue="protein_type", split=True, data=df, linewidth=0.2)
#[tick.set_rotation(90) for tick in ax.get_xticklabels()]
#plt.ylabel("Metabolic Function")
#plt.xlabel("C")
#plt.legend(title="Protein Type")
#plt.tight_layout()
#plt.show()
#plt.savefig("function_bdm_split_outliers.pdf")
#quit()


# ======== BDM histograms colored by interactions ==========

# need to get better interaction data!

# -----> just the viruses, but color by host interaction taxonomy
df = df.loc[df['protein_type'] == 'virus']
tax = 'phylum'
ax = sns.kdeplot(data=df, x="whole_bdm", hue=tax, log_scale=(True, True), fill=False,
                 legend=False, alpha=.5, linewidth=0.7, palette="hls")
ax.set(xlabel='BDM Value', ylabel='Density', title='Host Proteins ' + tax)
plt.savefig('Figures/virus_bdm_'+tax+'.pdf')
plt.clf()

# -----> just the viruses, but color by host interaction taxonomy
df = df.loc[df['protein_type'] == 'virus']
tax = 'phylum'
ax = sns.kdeplot(data=df, x="whole_bdm", hue=tax, log_scale=(True, False), multiple="fill", linewidth=0.1, legend=False,
                 palette="hls")
ax.set(xlabel='BDM Value', ylabel='Density', title='Host Proteins ' + tax)
plt.savefig('Figures/virus_bdm_'+tax+'_fraction.pdf')
plt.clf()
quit()


# ======== BDM vs degree =============

#sns.scatterplot(x='whole_bdm', y='degree', data=df, s=2, alpha=0.7, hue='protein_type', linewidth=0)  #, ncolors=nfam
#plt.xscale('log')
#plt.xlabel("C")
#plt.ylabel("degree")
#plt.tight_layout()
#plt.show()
#plt.savefig('bdm_degree.pdf')

# These are a small group of proteins near the top
#host_points = df.loc[(df['protein_type'] == 'host') & (df['degree'] > 400) & (df['degree'] < 550) & (df['whole_bdm'] > 1100) & (df['whole_bdm'] < 1400)]

# bdm vs degree host proteins only, colored by taxonomy TODO: Fix the degree :(
#df = df.loc[df.degree < 800]  # remove outliers
#df = df.loc[df['protein_type'] == 'host']
#df = df.dropna(subset=['kingdom'])  # get rid of nans for interested column
#df['kingdom'] = df['kingdom'].astype(int)  # make tax_ids integers
#ncolors = len(list(set(df['kingdom'])))  # get the number of distinct colors
#sns.set_palette("YlGnBu", n_colors=ncolors)
#sns.scatterplot(x='whole_bdm', y='degree', data=df, s=1, alpha=0.7, hue='kingdom', linewidth=0, palette="colorblind", legend="full")  #, ncolors=nfam
#plt.xscale('log')
#plt.xlabel("C")
#plt.ylabel("Degree")
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1, fontsize=8, title="Kingdom ID", title_fontsize=9)
#plt.tight_layout()
#plt.show()
#plt.savefig('bdm_degree_kingdom.pdf')

# bdm vs degree virus proteins only, colored by species  # TODO: FIX THE DEGREE :(
#df = df.loc[df['protein_type'] == 'virus']
#sns.scatterplot(x='whole_bdm', y='degree', data=df, s=1, alpha=0.7, hue='species', linewidth=0, palette="colorblind", legend="full")  #, ncolors=nfam
#plt.xscale('log')
#plt.xlabel("C")
#plt.ylabel("Degree")
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1, fontsize=8, title="Species", title_fontsize=9)
#plt.tight_layout()
#plt.show()
#plt.savefig('bdm_degree_virus_species.pdf')

# Avg BDM vs component size (nothing)
#comps = nx.connected_components(graph)
#sizes = []
#avg_bdms = []
#for comp in comps:
#    if len(comp) < 100:
#        continue
#    sizes.append(len(comp))
#    bdms = list(map(lambda x: graph.nodes[x]['whole_bdm'], list(comp)))
#    avg_bdms.append(sum(bdms)/len(bdms))
#plt.scatter(sizes, avg_bdms)
#plt.show()

# components
#n_comps = nx.number_connected_components(graph)
#comps = nx.connected_components(graph)
#nodes_in_comps = list(map(lambda x: len(x), list(comps)))
#nodes_in_comps.sort(reverse=True)
#plt.plot(nodes_in_comps)
#plt.xscale("log")
#plt.yscale("log")
#plt.savefig("component_node_dist.pdf")
#quit()

# largest connected component-- Need to work on this
#largest_cc = max(comps, key=len)
#largest_cc = graph.subgraph(largest_cc).copy()
#nx.draw(largest_cc)
#plt.savefig("lcc.pdf")
#quit()

# Degree distribution histogram
#hist = nx.degree_histogram(graph)
#plt.plot(hist)
#plt.xscale("log")
#plt.yscale("log")
#plt.savefig("degree_dist.pdf")


# ================ Try to visualize the PPI network ===================

#plt.figure(1, figsize=(8, 8))
# layout graphs with positions using graphviz neato
#pos = graphviz_layout(graph, prog="neato")
# color nodes the same in each connected subgraph
#C = (graph.subgraph(c) for c in nx.connected_components(graph))
#for g in C:
#    print(g)
#    c = [random.random()] * nx.number_of_nodes(g)  # random color...
#    nx.draw(g, pos, node_size=40, node_color=c, vmin=0.0, vmax=1.0, with_labels=False)
#plt.show()

#figure(figsize=(10, 8))
#pos = nx.spring_layout(graph, iterations=10)  # with_labels=True
#nx.draw(graph, pos)
#plt.show()




'''
print('step 2')

# Load in the PPIs
# With each line, get the different in BDM
first_time_running = True

bdm_diffs = []

if first_time_running:
    with open(os.path.join(ppis, 'protein.links.full.v10.5.txt'), 'rb') as f:

        # skip first line, save header
        header = next(f)
        header = header.split()
        n = 0

        # read in line by line
        for i, line in enumerate(f):

            # reduce the dataset TODO: Implement dask!
            # random 1 in 20 chance it gets included in plot
            #if random.choice(range(20)) != 0:
            #    continue

            line = line.split()

            # only load in the experiment and database ones
            experiments = int(line[9])
            database = int(line[11])
            if experiments == 0 and database == 0:
                continue

            p1 = line[0].decode("utf-8")
            p2 = line[1].decode("utf-8")

            # might not actually be in the bdm dict, if not then continue
            try:
                p1bdm = bdm_dict[p1]['whole_bdm']
                p2bdm = bdm_dict[p2]['whole_bdm']
                diff = p1bdm - p2bdm
                bdm_diffs.append(diff)
                n += 1
            except:
                continue

            # pickle all the bdm diffs so we don't recalculate
            if n % 1000000 == 0:
                with open(os.path.join(pickle_out, str(int(n/1000000)) + '_bdm_diffs.p'), 'wb') as handle:
                    pickle.dump(bdm_diffs, handle)
                    bdm_diffs = []
                    print(n)


# load in all the BDM differences
files = os.listdir(pickle_out)
files = list(filter(lambda x: not re.search('DS_Store', x) and not re.search('rand', x) and not re.search('whole', x), files))

bdm_diffs = []

for i, file in enumerate(files):

    # too much to plot! Runs out of application memory. Halve the data?
    with open(os.path.join(pickle_out, file), 'rb') as f:
        bdms = pickle.load(f)
        bdm_diffs.extend(bdms)
        print(file)

bdm_diffs = [abs(x) for x in bdm_diffs]
bdm_diffs.sort(reverse=True)


# ==================================================================================================
# --------------------------- Step 3: Calculate BDM diffs from randomly made links ---------------------------
# ==================================================================================================

print('step 3')

# There are 600M pairs
# About 10M proteins
# need to make a randomized sample, several times, for comparison
# Can also solve analytically

n_dists = 100
n_pairs = len(bdm_diffs)

first_time_running = False

if first_time_running:
    for dist in range(n_dists):

        # select random proteins
        random_pairs = [random.choice(bdm_dist) for _ in range(n_pairs*2)]

        # make them into pairs and get their bdm diffs
        random_bdm_diffs = zip(*[iter(random_pairs)]*2)
        random_bdm_diffs = map(lambda x: abs(x[0] - x[1]), random_bdm_diffs)
        random_bdm_diffs = list(random_bdm_diffs)
        random_bdm_diffs.sort(reverse=True)

        # pickle them for later because of memory
        with open(os.path.join(pickle_out, str(dist) + '_rand_bdm_diffs.p'), 'wb') as handle:
            pickle.dump(random_bdm_diffs, handle)
        print(dist)


# load in all the random-pair BDM differences
files = os.listdir(pickle_out)
files = list(filter(lambda x: not re.search('DS_Store', x) and re.search('rand', x) and not re.search('whole', x), files))[:5]

rand_bdm_diffs = []

for i, file in enumerate(files):
    # too much to plot! Runs out of application memory. Halve the data?
    with open(os.path.join(pickle_out, file), 'rb') as f:
        bdms = pickle.load(f)
        bdms = [abs(x) for x in bdms]
        bdms.sort(reverse=True)
        rand_bdm_diffs.append(bdms)
        print(file)


# plot them all on a single plot
for d in rand_bdm_diffs:
    plt.plot(d, color='grey', alpha=0.3)
plt.plot(bdm_diffs, color='red', alpha=1)
#plt.yscale('log')
plt.ylabel(r'$\Delta$ BDM')
plt.savefig('random_ppis.png')
plt.show()
'''