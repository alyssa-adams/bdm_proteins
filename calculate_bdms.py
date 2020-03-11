import os
from functions import Complexity


# this script calculates all the bdm values for all possible window sizes for

#type = 'proteins'
type = 'random_proteins'
grouping = '9'
pickle_out = 'bdms_' + type + '_' + grouping
pickle_out = os.path.join('pickle_jar', pickle_out)

# initialize classes
genecomplexity = Complexity()
bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

genecomplexity.calculate_bdms(bdm=bdm, pickle_out=pickle_out, type=type, grouping=grouping, data_directory='t7')
