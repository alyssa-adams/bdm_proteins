import os
from functions import Complexity


'''
This script calculates all the bdm values for window sizes 3-12

Accepted types: proteins, random_proteins, genes
Accepted groupings: EDSSMat90, 9
'''

type = 'genes'
grouping = ''
data_directory = 'data_covid_evolution'
pickle_out = 'bdms_' + type + '_' + grouping + '_' + data_directory
pickle_out = os.path.join('bdm_pickles', pickle_out)

# initialize classes
genecomplexity = Complexity()
bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

# calculate it and save to the pickle file
genecomplexity.calculate_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory=data_directory)
