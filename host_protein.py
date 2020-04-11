from functions import Complexity
import os, re, csv, pickle
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')




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
