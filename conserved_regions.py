import os, csv
from functions import Complexity
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt


'''
File:
autographivirinae_tail-fiber_by-host.motif-seqs (1).tsv

End of the tail:
QVWSGSAGGGVSVTVSQDLRFRNIWIKCANNSWNFFRTGPDGIYFIASDGGWLRFQIHSNGLGFKNIADSRSVPNAIMVENE

Get all the conserved regions and their variances
dict = {region: [list of variances]}

Measure BDM of the base region, then get the BDM for the variants

Also get the BDM for all the other regions along that tail that are NOT conserved

2 plots: 
Violin plot of variant regions and their BDM values with one red dot indicating the BDM of the base region.
Two violins, one of all the base regions and the other of the regions that are not conserved.
'''


tail = 'QVWSGSAGGGVSVTVSQDLRFRNIWIKCANNSWNFFRTGPDGIYFIASDGGWLRFQIHSNGLGFKNIADSRSVPNAIMVENE'

# Let's use EDSSMat90 grouping
grouping = 'EDSSMat90'
type = 'proteins'

# initialize classes
genecomplexity = Complexity()
bdm = genecomplexity.init_bdm(type=type, grouping=grouping)

# read in file
in_filepath = os.path.join('conserved_regions_data', 'autographivirinae_tail-fiber_by-host.motif-seqs (1).tsv')

conserved_regions = {}

with open(in_filepath) as csvfile:

    reader = csv.reader(csvfile, delimiter='\t')

    # skip header
    header = next(reader)

    for row in reader:

        base_region = row[-2]
        variant_region = row[-1]

        if base_region in conserved_regions.keys():
            conserved_regions[base_region].append(variant_region)
        else:
            conserved_regions[base_region] = [variant_region]


# now get all the 6-length regions that are NOT variable
not_conserved_regions = []

for i, letter in enumerate(tail[:-5]):

    region = tail[i:i+6]

    if region not in conserved_regions.keys():
        not_conserved_regions.append(region)


# --------------- FOR PLOT ONE ---------------

#Measure BDM of the base region, then get the BDM for the variants
# sequence: base_bdm, average_variable_bdm, bdm for each variant

conserved_bdms = {}

for sequence in conserved_regions.keys():

    conserved_bdms[sequence] = {}

    # base conserved bdm region
    base_seq = genecomplexity.group_amino_acids(grouping, sequence)
    measure = bdm.bdm(base_seq, normalized=True)
    conserved_bdms[sequence]['base_bdm'] = measure

    # variants
    variants = conserved_regions[sequence]
    # convert sequences to 9 symbols
    variants = [genecomplexity.group_amino_acids(grouping, v) for v in variants]
    # get bdm for each
    variant_bdms = [bdm.bdm(v, normalized=True) for v in variants]
    conserved_bdms[sequence]['variant_bdms'] = dict(zip(conserved_regions[sequence], variant_bdms))

    # average of variants
    avg = sum(variant_bdms)/len(variant_bdms)
    conserved_bdms[sequence]['variant_bdms_avg'] = avg


# --------------- FOR PLOT TWO ---------------

#Also get the BDM for all the other regions along that tail that are NOT conserved
not_conserved_bdms = {}

for sequence in not_conserved_regions:

    base_seq = genecomplexity.group_amino_acids(grouping, sequence)
    measure = bdm.bdm(base_seq, normalized=True)
    not_conserved_bdms[sequence] = measure


# --------------- PLOT ONE ---------------

#Violin plot of variant regions and their BDM values with one red dot indicating the BDM of the base region.

x = [[key] * len(conserved_bdms[key]['variant_bdms']) for key in conserved_bdms.keys()]
x = [item for sublist in x for item in sublist]
y = [conserved_bdms[key]['variant_bdms'] for key in conserved_bdms.keys()]
y = [list(y[i].values()) for i in range(len(y))]
y = [item for sublist in y for item in sublist]
sns.violinplot(x=x, y=y)

plt.ylabel('BDM')
#plt.figure(figsize=(20, 4))
plt.xticks(rotation=65, horizontalalignment='right')
plt.show()

# folder stuff
figure_folder = 'conserved_regions_figures'
if not os.path.exists(figure_folder):
    os.makedirs(figure_folder)

plt.savefig(os.path.join(figure_folder, 'conserved_variants_bdm.png'))


# --------------- PLOT TWO ---------------

#Two violins, one of all the base regions and the other of the regions that are not conserved.

x = ['Conserved'] * len(conserved_bdms.keys()) + ['Not Conserved'] * len(not_conserved_bdms.keys())
y = [conserved_bdms[s]['base_bdm'] for s in conserved_bdms.keys()] + list(not_conserved_bdms.values())
sns.violinplot(x=x, y=y)
plt.ylabel('BDM')

# folder stuff
figure_folder = 'conserved_regions_figures'
if not os.path.exists(figure_folder):
    os.makedirs(figure_folder)

plt.savefig(os.path.join(figure_folder, 'conserved_notconserved_bdm.png'))


# variable regions of tails:
# variable bdm values, compare with random ones of the same size (update figure)
# Do variable regions have higher bdm value than non-variable ones?
# Do regions that are highly conserved have higher bdm values?
# What is the relationship between number of edits and bdm value?

# Plot: Number of edits (lev distance) and difference in BDM
