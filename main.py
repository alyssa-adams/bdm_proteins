# TODO: Need to make sure all the directories self-create
# TODO: Need to make sure all dependencies can auto-install, cut out any unnecessary ones
# TODO: Make the main script, so it can run from the command line

# assume data is downloaded in data folder
# assume PPIs are in data folder too
# For each protein, calculate the BDM (just do whole BDM for now)

from itertools import islice
import numpy as np
import os
import pickle, re
import pandas as pd
import random
from functions import Complexity
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('ggplot')
sns.set(rc={'figure.figsize': (11.7, 8.27)})


# specify file paths

type = 'proteins'
grouping = 'EDSSMat90'
data_source = 'virusstring'
ppis = os.path.join('data', data_source + '_ppis')
protein_sequences = os.path.join('data/', data_source + '_sequences')
pickle_out = os.path.join('pickle_jar/', data_source + '_bdms')


# ==================================================================================================
# --------------------------- Step 1: Get the BDMs of the proteins ---------------------------
# ==================================================================================================

print('step 1')

# initialize class from Functions.py file (only used to calculate bdms)
#genecomplexity = Complexity()
#bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

# run the calculation (takes hours, saves result to pickle files, only need to one once per set of proteins)
#genecomplexity.calculate_whole_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory=protein_sequences)

# or load in bdm values from pickle files from previous run
files = os.listdir(pickle_out)
files = list(filter(lambda x: not re.search('DS_Store', x) and not re.search('diffs', x), files))

# Save BDM values to df
bdm_dict = {}
for file in files:
    with open(os.path.join(pickle_out, file), 'rb') as f:
        bdms = pickle.load(f)
        bdm_dict.update(bdms)
df = pd.DataFrame.from_dict(bdm_dict, orient='index')
groups = list(set(df['group'].to_list()))
df = df.sort_values(by=['whole_bdm'])
df = df.reset_index(drop=True)
df = df.reset_index()

# Plot the distribution of BDMs
bdm_dist = list(map(lambda x: x[-1]['whole_bdm'], bdm_dict.items()))
bdm_dist.sort(reverse=True)
#ax = sns.scatterplot(x="index", y="whole_bdm", hue="group", data=df, linewidth=0, alpha=0.2, s=5, legend=False)
#plt.yscale('log')
#plt.gray()
#plt.ylabel('BDM')
#plt.savefig('group.pdf')  # , format='eps'
#plt.show()


# ==================================================================================================
# --------------------------- Step 2: Calculate BDM diffs from links ---------------------------
# ==================================================================================================

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
