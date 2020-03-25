BDM (Block Decomposition Method) is an approximation for algorithmic complexity. To see how this measure works, please see this nice demo and explanation: https://complexitycalculator.com/index.html

The BDM values are calculated using the PyBDM package: https://github.com/sztal/pybdm
The documentation for PyBDM is here: https://pybdm-docs.readthedocs.io/en/latest/#
And you can install PyBDM using conda.

This code:
1. Measures the complexity of biological sequences (proteins and genomes), and 
2. characterizes the complexity along the length of the sequence.

Scripts:

calculate_bdms.py - Run this to do the main calculations.
Input: folder full of .faa or .fna files
Output: pickle file of BDM measures (saved as a dict) in a separate directory bdm_pickles

functions.py - Contains the Complexity class, contains all the functions to process these files and do the calculations

main.py - Will be built to automate scripts they this can be run using terminal commands

sequence_heatmaps.py - Make one heatmap of complexity along the length of a sequence, one figure per sequence.
Input: bdm pickle file
Output: One figure per sequence in figures_whole_proteins

whole_bdms_figure.py - Makes a violin plot of BDM complexity measures, one measure per sequence, grouped by something.
Input: bdm pickle file
Output: One figure in whole_bdms_figure.py

conserved_regions.py - Script that explores conserved regions of a tail that is tracked with mutations
Input: a .tsv file of conserved regions and their variants
Output: violin plots


This repo is undergoing rapid development and will change a lot! 

The goal of this repo is this:
Given a set of proteins or genomes, see if BDM can provide insight into how proteins interact. 


