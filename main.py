# TODO: Need to make sure all the directories self-create
# TODO: Need to make sure all dependencies can auto-install, cut out any unnecessary ones
# TODO: Make the main script, so it can run from the command line

# assume data is downloaded in data folder
# assume PPIs are in data folder too
# For each protein, calculate the BDM (just do whole BDM for now)

import os, pickle, re
import random
from functions import Complexity
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt


type = 'proteins'
grouping = 'EDSSMat90'
data_source = 'virusstring'
ppis = os.path.join('data', data_source + '_ppis')
protein_sequences = os.path.join('data', data_source + '_sequences')
pickle_out = os.path.join('pickle_jar', data_source + '_bdms')

# initialize classes
#genecomplexity = Complexity()
#bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

# calculate it and save to the pickle file
#genecomplexity.calculate_whole_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory=protein_sequences)


# specify file stuff, get list of relevant files
files = os.listdir(pickle_out)
files = list(filter(lambda x: not re.search('DS_Store', x), files))[0:1]

'''
# Load in all the BDMs
bdm_dict = {}
for file in files:
    with open(os.path.join(pickle_out, file), 'rb') as f:
        bdms = pickle.load(f)
        bdm_dict.update(bdms)


# Load in the PPIs
# With each line, get the different in BDM

bdm_diffs = []

with open(os.path.join(ppis, os.listdir(ppis)[0]), 'rb') as f:

    # skip first line
    next(f)
    n = 0

    # all of this takes up way too much disk space
    # the PPI file alone is 31GB and has 600M+ rows
    # Instead, take a random sample of a 20th of these to plot

    # read in line by line
    for i, line in enumerate(f):

        # random 1 in 20 chance it gets included in plot
        if random.choice(range(20)) != 0:
            continue

        p1 = line.split()[0].decode("utf-8")
        p2 = line.split()[1].decode("utf-8")

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
'''


# load in all the BDM differences

files = os.listdir(pickle_out)
files = list(filter(lambda x: not re.search('DS_Store', x) and re.search('bdm_diffs', x), files))

bdm_diffs = []

for i, file in enumerate(files):
    # too much to plot! Runs out of application memory. Halve the data?
        with open(os.path.join(pickle_out, file), 'rb') as f:
            bdms = pickle.load(f)
            bdm_diffs.extend(bdms)

bdm_diffs = [abs(x) for x in bdm_diffs]
bdm_diffs.sort()

plt.plot(bdm_diffs)
plt.ylabel('BDM(Protein 1) - BDM(Protein 2)')
plt.savefig('virusstring_ppis.png')
plt.show()
