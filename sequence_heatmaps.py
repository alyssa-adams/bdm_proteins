from functions import Complexity
import os, re, csv, pickle
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.style.use('ggplot')
import numpy as np
import seaborn as sns; sns.set()


'''
This script makes a BDM heatmap figure per protein

Accepted types: proteins, random_proteins, genes
Accepted groupings: EDSSMat90, 9
'''

"""
# DEPRECIATED
def group_lineplots(bdm_file):

    '''
    Make a bunch of plots, one per tail
    :param bdm_file: str, the pickle file
    :param type: str, either genes or proteins
    :param grouping: str, 9, edsm9, 8
    :return:
    '''

    # folder stuff
    figure_folder = os.path.join('figures_old2', group)
    if not os.path.exists(figure_folder):
        os.makedirs(figure_folder)

    # load in the bdm values
    bdms_dict = pickle.load(open(bdm_file, 'rb'))
    type = bdm_file.split('_')[1]
    grouping = bdm_file.split('_')[2]
    window_size = bdm_file.split('_')[3]

    # restructure dict by groups
    bdms_dict_groups = {}

    for tail in bdms_dict.keys():
        group = bdms_dict[tail]['group']
        whole_bdm = bdms_dict[tail]['whole_bdm']
        bdms = bdms_dict[tail]['bdms']
        if len(bdms) == 0:
            continue

        if group not in bdms_dict_groups.keys():
            bdms_dict_groups[group] = {}

        bdms_dict_groups[group][tail] = (whole_bdm, bdms)

    # one plot per group
    for group in bdms_dict_groups.keys():

        # this is where the figure will go
        figure_folder = os.path.join('figures_old2')
        figure_out = os.path.join(figure_folder, group + '_' + type + '_' + grouping + '_' + window_size + '.png')
        #figure_folder = os.path.join('figures_old2', bdms_dict[tail]['group'])
        #figure_out = os.path.join(figure_folder, tail + '_' + type + '_' + grouping + '_' + window_size + '.png')

        # relevant data_tails for the plot
        tails = bdms_dict_groups[group]

        # get the max x length to normalize all the other lines to
        max_x = max([len(tails[tail][1]) for tail in tails.keys()])
        x_s = {}
        for tail in tails.keys():
            x1 = len(tails[tail][1])
            spacer = (max_x - x1)/x1
            x2 = [x + x*spacer for x in range(x1)]
            x_s[tail] = x2

        # initialize plot
        fig = plt.figure()
        ax = plt.axes()
        title = group + '_' + type + '_' + grouping + '_' + window_size
        plt.title(title)

        # Color the lines so that the bigger window sizes get darker colors
        #color_idx = np.linspace(0, 1, len(window_sizes))
        # , color=plt.cm.binary(color_idx[i])

        # Plot BDM
        for tail in tails.keys():
            ax.plot(x_s[tail], tails[tail][1], linewidth=0.5, alpha=0.7)

        # get whole measure
        #ax.axhline(y=whole_bdm, color='r', linewidth=0.6)

        plt.xlabel('BDM')
        plt.xlabel('Window (' + window_size + ')')
        ax.set_aspect('auto')
        plt.savefig(figure_out)


def normalize_x_scales(array):

    '''
    You have a bunch of lines, but you need to stretch them out along the same x scale
    :param array: Array of initial BDM values
    :return: list of lists x_s, of all x points for each line
    '''

    # get the ACTUAL PROTEIN SIZE length to normalize all the other lines to
    max_x = max([len(row) for row in array])
    # TODO: This depends on the smallest window size being 3
    protein_len = max_x + 2
    x_s = []
    for row in array:
        x1 = len(row)
        spacer = (protein_len - x1) / x1
        x2 = [x + x * spacer for x in range(x1)]
        x_s.append(x2)

    return x_s


def stretched_array(array, x_s):

    '''
    You have an array, but you want to normalize it by stretching the values across the whole array
    :param array: Array of initial BDM values
    :param x_s: Normalized x values, list of lists
    :return: stretched array
    '''

    # TODO: This depends on the smallest window size being 3
    max_x = max([len(row) for row in array])
    protein_len = max_x + 2

    # these are the bins to fit to
    bins = range(protein_len)
    new_array = []

    for i, x_1 in enumerate(x_s):

        x_values1 = x_1
        y_values1 = array[i]
        y_values2 = []

        # round DOWN the x_values1 to the nearest int as initial binning
        x_values1 = list(map(lambda x: int(x), x_values1))
        # these get to act like when the y values will "switch"

        for i, x1 in enumerate(x_values1):

            # if its the last one, take care of the edge case
            if i == len(x_values1) - 1:
                y2_chunk = [y_values1[i]] * (len(bins) - x_values1[i])
            else:
                y2_chunk = [y_values1[i]] * (x_values1[i + 1] - x_values1[i])
            y_values2 = y_values2 + y2_chunk

        new_array.append(y_values2)

    return new_array
"""

def right_padded_array(array):

    '''
    You have a list of lists with different row sizes. Align by padding with values at the end of each row.
    :param array: Array of initial BDM values
    :return: Aligned and padded array
    '''

    new_array = []
    max_x = max([len(row) for row in array])

    # For each row, just pad with the rest of the values
    for row in array:
        row2 = row + [row[0]] * (max_x - len(row))
        new_array.append(row2)

    return new_array


def sequence_heatmaps(bdm_file, figure_dir):

    '''
    Make a bunch of plots, one per tail
    :param bdm_file: str, the pickle file
    :param type: str, either genes or proteins
    :param grouping: str, 9, edsm9, 8
    :return:
    '''

    # folder stuff
    figure_folder = os.path.join(figure_dir)
    if not os.path.exists(figure_folder):
        os.makedirs(figure_folder)

    # load in the bdm values
    bdms_dict = pickle.load(open(bdm_file, 'rb'))

    # for each protein, make a single plot
    # for each protein, make an array with different sizes
    for protein in bdms_dict.keys():

        # get grouping
        #grouping = bdms_dict[protein]['group']
        grouping = 'EDSSMat90'

        # here is the array
        array = list(bdms_dict[protein]['bdms'].values())
        array.sort(key=len, reverse=True)

        # save the figure to this file
        figure_out = os.path.join(figure_folder, str(protein) + '_' + str(grouping) + '.png')

        # get the ACTUAL PROTEIN SIZE length to normalize all the other lines to
        #x_s = normalize_x_scales

        # to make a rectangular array for a heatmap, need to fill in the other points to make all rows the same size
        #new_array = stretched_array(array, x_s)

        # align all the values in columns by padding rows with the last value on the right size
        new_array = right_padded_array(array)

        # make about 15 tick marks
        x_n_ticks = 20
        x_ticks_space = int(len(new_array[0])/x_n_ticks)

        # the index of the position of yticks
        # TODO: yticklabels depends on min window size
        y_n_ticks = 3
        y_ticks_space = int(len(new_array) / y_n_ticks)
        y_values = list(range(3, len(new_array) + 3, y_ticks_space))
        # the spacing between ticks
        yticks = list(range(0, len(new_array), y_ticks_space))

        # initialize plot
        sns.set(font_scale=2)
        plt.subplots(figsize=(20, 4))
        ax = sns.heatmap(new_array, xticklabels=x_ticks_space, yticklabels=y_values, square=False,
                         cbar_kws={"orientation": "horizontal", "label": "BDM"}, cmap="Blues")

        # ticks
        ax.tick_params(axis='both', which='both', length=2)
        ax.set_yticks(yticks)
        ax.xaxis.tick_top()  # x axis on top
        ax.xaxis.set_label_position('top')

        # column lines
        ax.vlines(range(0, len(new_array[0])+1, 5), *ax.get_ylim(), linewidth=0.5)

        plt.xlabel('BP')
        plt.ylabel('Window Size')
        plt.title(str(protein) + ', AA Grouping: ' + str(grouping))
        #plt.title(str(protein), size=75, pad=40)
        plt.tight_layout()
        #plt.show()
        plt.savefig(figure_out, format='png')


type = 'proteins'
grouping = 'EDSSMat90'

# input file
#bdm_file_name = 'bdms' + '_' + type + '_' + grouping
bdm_file_name = 'bdms_proteins_EDSSMat90_data_t7_host'
bdm_file = os.path.join('bdm_pickles', bdm_file_name)

# save these figures to these directories
main_dir = 'figures_whole_proteins'
type_dir = type
grouping_dir = os.path.join(grouping, 't7_host')
#grouping_dir = 'covid'
figure_dir = os.path.join(main_dir, type_dir, grouping_dir)

# Make the figures
sequence_heatmaps(bdm_file=bdm_file, figure_dir=figure_dir)
