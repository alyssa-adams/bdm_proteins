import re, os
from pybdm import BDM
from pybdm.partitions import PartitionRecursive
import numpy as np
import pickle
import random
from snapgene_reader import snapgene_file_to_dict, snapgene_file_to_seqrecord


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
        if type == 'proteins' or type == 'random_proteins':
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


    # make randomly generated proteins
    def make_random_proteins(self, max_size, number_of_proteins):

        '''
        Make some random proteins
        :param max_size: Max size of the proteins
        :param number_of_proteins: The number of them
        :return: Dict {id: group, len, sequence}
        '''

        proteins = {}

        for protein in range(number_of_proteins):

            # TODO: Min protein size of 80
            size = random.choice(range(80, max_size))
            sequence = [random.randint(0, 8) for _ in range(0, size)]

            proteins[protein] = {
                # TODO: Random groups?
                'group': random.choice(range(0, 10)),
                'length': len(sequence),
                'sequence': sequence
            }

        return proteins


    # translate amino acid sequences into 9 symbols
    def group_amino_acids(self, grouping, sequence):

        '''
        Group the amino acids together into 9 symbols
        :param grouping: str, either 9 or EDSSMat90
        :param sequence: str, the sequence to translate
        :return: np array, the translated sequence
        '''

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

            sequence = re.sub('[G]{1}', '0', sequence)
            sequence = re.sub('[P]{1}', '1', sequence)
            sequence = re.sub('[ALVI]{1}', '2', sequence)
            sequence = re.sub('[CM]{1}', '3', sequence)
            sequence = re.sub('[DE]{1}', '4', sequence)
            sequence = re.sub('[RK]{1}', '5', sequence)
            sequence = re.sub('[FYWH]{1}', '6', sequence)
            sequence = re.sub('[NQ]{1}', '7', sequence)
            sequence = re.sub('[TS]{1}', '8', sequence)

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

            sequence = re.sub('[ILVATM]{1}', '0', sequence)
            sequence = re.sub('[RK]{1}', '1', sequence)
            sequence = re.sub('[N]{1}', '2', sequence)
            sequence = re.sub('[DE]{1}', '3', sequence)
            sequence = re.sub('[CFWY]{1}', '4', sequence)
            sequence = re.sub('[P]{1}', '5', sequence)
            sequence = re.sub('[S]{1}', '6', sequence)
            sequence = re.sub('[G]{1}', '7', sequence)
            sequence = re.sub('[HQ]{1}', '8', sequence)

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

            sequence = re.sub('[AST]{1}', '0', sequence)
            sequence = re.sub('[RQKHE]{1}', '1', sequence)
            sequence = re.sub('[ND]{1}', '2', sequence)
            sequence = re.sub('[C]{1}', '3', sequence)
            sequence = re.sub('[ILMV]{1}', '4', sequence)
            sequence = re.sub('[FYW]{1}', '5', sequence)
            sequence = re.sub('[P]{1}', '6', sequence)
            sequence = re.sub('[G]{1}', '7', sequence)

        else:
            print('No amino acid grouping specified!')

        sequence = np.array(list(sequence), dtype=int)

        return sequence


    def translate_sequence(self, sequence, grouping):

        '''
        translate the sequence of base pairs or amino acids into a np array
        :param sequence: string, the sequence to translate
        :param grouping: str, proteins or nuc
        :return: sequence of ints, np array
        '''

        sequence = sequence.upper()

        if type == 'proteins':

            if re.search('>', sequence):
                sequence
            sequence = self.group_amino_acids(grouping, sequence)

        else:

            sequence = re.sub('A', '0', sequence)
            sequence = re.sub('T', '1', sequence)
            sequence = re.sub('C', '2', sequence)
            sequence = re.sub('G', '3', sequence)

            sequence = np.array(list(sequence), dtype=int)

        return sequence


    # groupings: 9, EDSSMat90, 8
    def read_file(self, path, type, grouping):

        '''
        Read in genome file and save to dict
        :param path: File path
        type: str, either 'genes' or 'proteins'
        grouping: str, the particular amino acid grouping. None if not using proteins.
        :return: Dict {id: len, sequence}
        '''

        sequences = {}


        # check to see if this is a .dna file type

        extension = os.path.splitext(path)[1]
        id = os.path.basename(path).split('.')[0]

        if extension == '.dna':
            file_dict = snapgene_file_to_dict(path)
            sequence = file_dict['seq']

            # translate into ints
            sequence = self.translate_sequence(sequence, grouping)

            sequences[id] = {
                'group': id,
                'length': len(sequence),
                'sequence': sequence
            }

            return sequences


        # read in the whole file
        with open(path) as file:
            content = file.read()

        # get the proteins in text blobs using regex
        whole_regex = '>[^\>]*'
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
            # TODO: Check threshold of aa sequences that are too small. Now it's at 80.
            if len(sequence) < 80:
                continue

            # throw away sequences with an X in them
            if re.search('X', sequence):
                continue

            # N = Any of the ATCG
            # Trim away anything that isn't ATCG
            if re.search('[^ATCG]', sequence):
                sequence = re.sub('[^ATCG]', '', sequence)

            # translate into ints
            sequence = self.translate_sequence(sequence, grouping)

            sequences[id] = {
                'group': group,
                'length': len(sequence),
                'sequence': sequence
            }

        return sequences


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
        # TODO: The smallest window size is currently 3aa, DON'T go past 12. Just to up through 12.
        window_sizes = range(3, 6)  #TODO Change this back!!!! And the things below too!!!!
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


    def all_bdms(self, bdm, sequence, grouping):

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

        # integrate BDM
        bdm_list = [item for sublist in list(bdms.values()) for item in sublist]
        bdms_mass = sum(bdm_list)
        bdms_area = len(bdm_list)
        bdm_density = bdms_mass / bdms_area

        '''
        # --- second-order bdm ---

        # for each key, round the values to the grouping
        if grouping == 'EDSSMat90' or grouping == '9':
            n_bins = 9
        else:
            n_bins = 4

        # make the bins
        bin_size = 1/float(n_bins)
        bins = tuple(zip([i*bin_size for i in range(n_bins)], [i*bin_size + bin_size for i in range(n_bins)]))
        symbols = range(n_bins)
        bin_values = [i * bin_size + bin_size for i in range(n_bins)]
        dictionary = dict(zip(bin_values, symbols))

        # bin the bdm values
        second_order = []

        for key in bdms.keys():
            bdms_binned = []

            for value in bdms[key]:
                for bin in bins:
                    if value >= bin[0] and value < bin[1]:
                        bdms_binned.append(bin[1])
                        break

            # turn values into symbols
            second_bdm = [dictionary[i] for i in bdms_binned]

            # the actual bdm value
            second_bdm = bdm.bdm(np.array(second_bdm), normalized=True)
            second_order.append(second_bdm)


        # --- third-order bdm ---

        # bin these values
        third_order = []
        for value in second_order:
            for bin in bins:
                if value >= bin[0] and value < bin[1]:
                    third_order.append(bin[1])
                    break

        # turned binned values into integers
        third_order = [dictionary[i] for i in third_order]
        third_order = bdm.bdm(np.array(third_order), normalized=True)

        # average second order?
        second_order = sum(second_order) / len(second_order)
        '''
        second_order = None
        third_order = None

        return {'whole_bdm': whole_bdm,
                'bdm_density': bdm_density,
                'second_order': second_order,
                'third_order': third_order,
                'bdms': bdms}


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

        # going to save the output in a dict format
        bdms_dict = {}

        # if random_proteins, then just make the random proteins instead of looping over files
        if type == 'random_proteins':

            proteins = self.make_random_proteins(max_size=200, number_of_proteins=100)

            # say how many proteins there are
            print(len(proteins.keys()))

            # for each protein, calculate all the bdms for all the windows
            for protein in proteins.keys():

                # say which protein its on
                print(protein)

                sequence = np.array(proteins[protein]['sequence'])

                # get the values for each frame, for each window size
                # bdms is a dict, where key is window size and value is list of values
                bdms = self.all_bdms(bdm, sequence, grouping=grouping)

                # save the values to the pickle dict
                bdms_dict[protein] = {
                    'group': proteins[protein]['group'],
                    'whole_bdm': bdms['whole_bdm'],
                    'bdm_density': bdms['bdm_density'],
                    'second_order': bdms['second_order'],
                    'third_order': bdms['third_order'],
                    'bdms': bdms['bdms']
                }

            # pickle to save the bdm values
            with open(pickle_out, 'wb') as handle:
                pickle.dump(bdms_dict, handle)

            return None


        # Else, read in the files

        # specify file stuff
        files = os.listdir(data_directory)
        files = list(filter(lambda x: not re.search('DS_Store', x), files))

        # check to see which files to look at
        if type == 'proteins':
            files = list(filter(lambda x: re.search('\.faa', x) or re.search('\.dna', x), files))
            # don't set a hard window_size here, so it can change

        elif type == 'genes':
            files = list(filter(lambda x: re.search('\.ffn', x) or re.search('\.dna', x) or re.search('\.fasta', x), files))

        else:
            print('Please specify type as either genes or proteins!')


        # for each file found in the directory
        for file in files:

            print(file)

            # specify paths
            sequence_file = os.path.join(data_directory, file)

            # read in the file
            sequences = self.read_file(sequence_file, type=type, grouping=grouping)

            # say how many proteins there are
            print(len(sequences.keys()))

            # for each protein, calculate all the bdms for all the windows
            for p in sequences.keys():

                # say which protein its on
                print(p)

                try:  # TODO: This needs to be fixed

                    sequence = sequences[p]['sequence']
                    group = sequences[p]['group']

                    # get the values for each frame, for each window size
                    # bdms is a dict, where key is window size and value is list of values
                    bdms = self.all_bdms(bdm, sequence, grouping=grouping)

                    # save the values to the pickle dict
                    bdms_dict = {
                        'group': group,
                        'whole_bdm': bdms['whole_bdm'],
                        'bdm_density': bdms['bdm_density'],
                        'second_order': bdms['second_order'],
                        'third_order': bdms['third_order'],
                        'bdms': bdms['bdms']
                    }

                    # pickle to save the bdm values
                    p = re.sub('/', '|', p)
                    folder = os.path.join(pickle_out)
                    if not os.path.exists(folder):
                        os.makedirs(folder)
                    with open(pickle_out + '/' + p, 'wb') as handle:
                        pickle.dump(bdms_dict, handle)

                except:
                    print('failed')
                    continue

        return None
