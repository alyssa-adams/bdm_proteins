# TODO: Need to make sure all the directories self-create
# TODO: Need to make sure all dependencies can auto-install, cut out any unnecessary ones
# TODO: Make the main script, so it can run from the command line

# assume data is downloaded in data folder
# assume PPIs are in data folder too
# For each protein, calculate the BDM (just do whole BDM for now)

import os, pickle
from functions import Complexity


type = 'proteins'
grouping = 'EDSSMat90'
data_source = 'virusstring'
ppis = os.path.join('data', data_source + '_ppis')
protein_sequences = os.path.join('data', data_source + '_sequences')
pickle_out = os.path.join('pickle_jar', data_source + '_bdms')

# initialize classes
genecomplexity = Complexity()
bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

# calculate it and save to the pickle file
genecomplexity.calculate_whole_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory=protein_sequences)

# load in pickle file
#bdms = pickle.load(open(os.path.join(pickle_out, 'whole_bdms.p'), 'rb'))

# Use the PPI file to plot!
#with open(os.path.join(ppis, 'protein.links.v10.5.txt'), 'r') as f:
#    for line in f:
#        p1 = line.split()[0]
#        p2 = line.split()[1]

#ppi_data
