# TODO: Need to make sure all the directories self-create
# TODO: Need to make sure all dependencies can auto-install, cut out any unnecessary ones
# TODO: Make the main script, so it can run from the command line

# assume data is downloaded in data folder
# assume PPIs are in data folder too
# For each protein, calculate the BDM (just do whole BDM for now)

import os
from functions import Complexity


type = 'proteins'
grouping = 'EDSSMat90'
ppis = 'data/virusstring_ppis'
protein_sequences = 'data/virusstring_sequences'
pickle_out = 'bdms_' + type + '_' + grouping + '_' + data_directory
pickle_out = os.path.join('bdm_pickles', pickle_out)

# initialize classes
genecomplexity = Complexity()
bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

# calculate it and save to the pickle file
genecomplexity.calculate_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory=data_directory)
