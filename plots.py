from functions import Complexity
import os, re, csv, pickle
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.style.use('ggplot')
import numpy as np


def plot_group(bdm_file):

    '''
    Make a bunch of plots, one per tail
    :param bdm_file: str, the pickle file
    :param type: str, either genes or proteins
    :param grouping: str, 9, edsm9, 8
    :return:
    '''

    # folder stuff
    figure_folder = os.path.join('figures2', group)
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
        figure_folder = os.path.join('figures2')
        figure_out = os.path.join(figure_folder, group + '_' + type + '_' + grouping + '_' + window_size + '.png')
        #figure_folder = os.path.join('figures2', bdms_dict[tail]['group'])
        #figure_out = os.path.join(figure_folder, tail + '_' + type + '_' + grouping + '_' + window_size + '.png')

        # relevant tails for the plot
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


def all_windows_3d(bdm_file):

    '''
    Make a bunch of plots, one per tail
    :param bdm_file: str, the pickle file
    :param type: str, either genes or proteins
    :param grouping: str, 9, edsm9, 8
    :return:
    '''

    # folder stuff
    figure_folder = os.path.join('figures2')
    if not os.path.exists(figure_folder):
        os.makedirs(figure_folder)

    # load in the bdm values
    bdms_dict = pickle.load(open(bdm_file, 'rb'))

    # for each protein, make a single plot
    # for each protein, make an array with different sizes
    for protein in bdms_dict.keys():

        # here is the array
        array = list(bdms_dict[protein]['bdms'].values())
        array.sort(key=len)

        # save the figure to this file
        figure_out = os.path.join(figure_folder, protein + '.png')

        # get the max x length to normalize all the other lines to
        max_x = max([len(row) for row in array])
        x_s = []
        for row in array:
            x1 = len(row)
            spacer = (max_x - x1)/x1
            x2 = [x + x*spacer for x in range(x1)]
            x_s.append(x2)

        # to make a rectangular array for a heatmap, need to fill in the other points to make all rows the same size
        for i in x_s:
            i
            # this is a test to make sure github is working

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


bdm_file = os.path.join('pickle_jar', 'bdms_proteins_EDSSMat90')
all_windows_3d(bdm_file=bdm_file)
