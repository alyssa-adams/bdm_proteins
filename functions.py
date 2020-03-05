import re, os
from pybdm import BDM
from pybdm.partitions import PartitionRecursive
import numpy as np
import pickle


# TODO: set up to run on many different kinds of proteins, known and unknown, for phages and hosts


class Complexity():

    '''
    This class reads in sequences and measures their complexity with BDM
    '''

    def init_bdm(self, type, grouping):

        '''
        Initialize the bdm class for genomes, always going to be ndim=1 and nsymbols=4
        type: str, either genes or proteins
        :return: bdm class
        '''

        # initialize BDM class
        if type == 'proteins':
            # only one amino acid grouping uses 8 symbols
            if grouping == '8':
                nsymbols = 8
            else:
                nsymbols = 9
        elif type == 'genes':
            nsymbols = 4
        else:
            print('What kind of sequence is it?')

        bdm = BDM(ndim=1, nsymbols=nsymbols, warn_if_missing_ctm=False, partition=PartitionRecursive)

        return bdm


    # groupings: 9, EDSSMat90, 8
    def read_file(self, path, type, grouping):

        '''
        Read in genome file and save to dict
        :param path: File path
        type: str, either 'genes' or 'proteins'
        grouping: str, the particular amino acid grouping. None if not using proteins.
        :return: Dict {id: len, sequence}
        '''

        proteins = {}

        # read in the whole file
        with open(path) as file:
            content = file.read()

        # get the proteins in text blobs using regex
        whole_regex = '>[^\*]*'
        parse_us = re.findall(whole_regex, content)

        for blob in parse_us:

            # if empty row, skip
            if len(blob) == 0:
                continue

            # split by newline, join all except the first
            blob_parts = blob.split('\n')

            # get id and group, save to dict
            id = re.sub('>', '', blob_parts[0].split()[0])
            group = blob_parts[0].split()[-1]

            # make sequence from the other parts
            sequence = ''.join(blob_parts[1:])

            # throw away sequences that are too small
            # TODO: Check threshold of aa sequences that are too small. Now it's at 10.
            if len(sequence) < 10:
                continue

            # throw away sequences with an X in them
            if re.search('X', sequence):
                continue

            if type == 'proteins':

                sequence = re.sub('\*', '', sequence)

                if grouping == '9':

                    # group the amino acids down to 9 groups
                    # G
                    # P
                    # ALVI
                    # CM
                    # DE
                    # RK
                    # FYWH
                    # NQ
                    # TS

                    sequence = re.sub('[G]{1,}', '0', sequence)
                    sequence = re.sub('[P]{1,}', '1', sequence)
                    sequence = re.sub('[ALVI]{1,}', '2', sequence)
                    sequence = re.sub('[CM]{1,}', '3', sequence)
                    sequence = re.sub('[DE]{1,}', '4', sequence)
                    sequence = re.sub('[RK]{1,}', '5', sequence)
                    sequence = re.sub('[FYWH]{1,}', '6', sequence)
                    sequence = re.sub('[NQ]{1,}', '7', sequence)
                    sequence = re.sub('[TS]{1,}', '8', sequence)

                elif grouping == 'EDSSMat90':

                    # put second grouping here!
                    # ILVATM
                    # RK
                    # N
                    # DE
                    # CFWY
                    # P
                    # S
                    # G
                    # HQ

                    sequence = re.sub('[ILVATM]{1,}', '0', sequence)
                    sequence = re.sub('[RK]{1,}', '1', sequence)
                    sequence = re.sub('[N]{1,}', '2', sequence)
                    sequence = re.sub('[DE]{1,}', '3', sequence)
                    sequence = re.sub('[CFWY]{1,}', '4', sequence)
                    sequence = re.sub('[P]{1,}', '5', sequence)
                    sequence = re.sub('[S]{1,}', '6', sequence)
                    sequence = re.sub('[G]{1,}', '7', sequence)
                    sequence = re.sub('[HQ]{1,}', '8', sequence)

                elif grouping == '8':

                    # put third grouping here!
                    # AST
                    # RQKHE
                    # ND
                    # C
                    # ILMV
                    # FYW
                    # P
                    # G

                    sequence = re.sub('[AST]{1,}', '0', sequence)
                    sequence = re.sub('[RQKHE]{1,}', '1', sequence)
                    sequence = re.sub('[ND]{1,}', '2', sequence)
                    sequence = re.sub('[C]{1,}', '3', sequence)
                    sequence = re.sub('[ILMV]{1,}', '4', sequence)
                    sequence = re.sub('[FYW]{1,}', '5', sequence)
                    sequence = re.sub('[P]{1,}', '6', sequence)
                    sequence = re.sub('[G]{1,}', '7', sequence)

                else:
                    print('No amino acid grouping specified!')

            else:

                sequence = re.sub('A', '0', sequence)
                sequence = re.sub('T', '1', sequence)
                sequence = re.sub('C', '2', sequence)
                sequence = re.sub('G', '3', sequence)

            sequence = np.array(list(sequence), dtype=int)

            proteins[id] = {
                'group': group,
                'length': len(sequence),
                'sequence': sequence
            }

        return proteins


    def bdm_slidingwindow(self, bdm, sequence):

        '''
        Do every possible window size
        Give it a window size, and then calculate the bdm for each window
        bdm: The initialized bdm class
        :param sequence: ATCGetc
        :return: A dict of bdm values, keys are window size and values are a list of bdms, one value per window
        '''

        bdms = {}

        # Calculate window sizes based on length of the sequence
        # All possible window sizes, starting at 5aa

        # make window sizes
        window_sizes = range(5, len(sequence))
        for size in window_sizes:

            measures = []

            for window in range(len(sequence) - size + 1):

                # make sure slice size isn't a fragment, if it is, just skip
                slice = sequence[window:window+size]
                if len(slice) < size:
                    continue
                measure = bdm.bdm(slice, normalized=True)
                measures.append(measure)

            bdms[size] = measures

        return bdms


    def all_bdms(self, bdm, sequence):

        '''
        Gets all bdm values for each window size, and also the entire string
        :param bdm: initialized bdm class
        :param sequence: whole genome/gene sequence
        :return: bdm for entire sequence, list of bdm values for the sliding window, and the window size it picked
        '''

        # BDM for the whole sequence
        whole_bdm = bdm.bdm(sequence, normalized=True)

        # do sliding window of Xbp until end
        bdms = self.bdm_slidingwindow(bdm, sequence)

        return whole_bdm, bdms


    def calculate_bdms(self, bdm, pickle_out, type, grouping, data_directory):

        '''
        Calculates all the BDM values for files of sequences
        :param bdm: initialized bdm class
        :param pickle_out: str, The file to save these values to
        :param type: str, either 'genes' or 'proteins'
        :param grouping: str, amino acid groups, either '9' 'EDSSMat90' or '8'
        :param data_directory: specify which directory to find all the raw sequences
        :return: None, just makes file
        '''

        # specify file stuff
        files = os.listdir(data_directory)
        files = list(filter(lambda x: not re.search('DS_Store', x), files))

        # check to see which files to look at
        if type == 'proteins':
            files = list(filter(lambda x: re.search('\.faa', x), files))
            # don't set a hard window_size here, so it can change

        elif type == 'genes':
            files = list(filter(lambda x: re.search('\.ffn', x), files))

        else:
            print('Please specify type as either genes or proteins!')

        # going to save the output in a dict format
        bdms_dict = {}

        # for each file found in the directory
        for file in files:

            print(file)

            # specify paths
            group = file.split('_')[0]
            proteins_file = os.path.join(data_directory, file)

            # read in the file
            proteins = self.read_file(proteins_file, type=type, grouping=grouping)

            # say how many proteins there are
            print(len(proteins.keys()))

            # for each protein, calculate all the bdms for all the windows
            for protein in proteins.keys():

                # say which protein its on
                print(protein)

                sequence = proteins[protein]['sequence']

                # get the values for each frame, for each window size
                # bdms is a dict, where key is window size and value is list of values
                whole_bdm, bdms = self.all_bdms(bdm, sequence)

                # save the values to the pickle dict
                bdms_dict[protein] = {
                    'group': group,
                    'whole_bdm': whole_bdm,
                    'bdms': bdms
                }

        # pickle to save the bdm values
        with open(pickle_out, 'wb') as handle:
            pickle.dump(bdms_dict, handle)

        return None
