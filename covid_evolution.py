from functions import Complexity
import pandas as pd
import os, re, csv, pickle
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import difflib
from datetime import datetime
from collections import Counter


# Just look at the aligned sequences
# separate by region, sort by date
# Just sort by date

# overall, how often does it mutate?
# overall, how often does it mutate PER REGION?

# Is there a subsequence where it mutates the most?

def make_dates_countries():

    # read in the aligned sequences
    type = 'genes'
    grouping = ''
    data_directory = 'data_covid_evolution_aligned'
    path = os.path.join(data_directory, 'covid_evolution_aligned.pir')
    genecomplexity = Complexity()
    data = genecomplexity.read_file(path, type, grouping, aligned=True)

    # remake the dict to get the dates and countries
    data_dates_countries = {}

    for d in data.keys():

        sequence = data[d]['sequence']

        group = data[d]['group']
        date = re.findall('\d{4}\-\d{2}\-\d{2}', group)

        # if there isn't a date, then we don't care
        if len(date) == 0:
            continue

        date = date[0]

        country = group.split('/')[1]

        data_dates_countries[d] = {
            'date': date,
            'country': country,
            'sequence': sequence
        }

    # save this dict to a pickle file
    with open(os.path.join('bdm_pickles', 'covid_dates_countries_dict'), 'wb') as handle:
        pickle.dump(data_dates_countries, handle)


#make_dates_countries()

data_dates_countries = pickle.load(open(os.path.join('bdm_pickles', 'covid_dates_countries_dict'), 'rb'))

# make into a df
df = pd.DataFrame.from_dict(data_dates_countries)
df = df.transpose()

# sort by country
df_countries = {}
for country, df_country in df.groupby('country'):
    df_country = df_country.sort_values(by=['date'])
    df_countries[country] = df_country


def mutations(df):

    mutations = {}

    # the list of the sequences
    sequences = df['sequence'].tolist()
    bps = len(sequences[0])

    for bp in range(bps):

        # for each bp location, see which ones are different
        location_bps = list(map(lambda x: x[bp], sequences))
        bps_freq = dict(Counter(location_bps))

        # if this is only one symbol, skip
        if len(bps_freq.keys()) < 2:

            mutations[bp] = {
                'bps_freq': None,
                'bps_freq_norm': None,
                'mutation_locations': None
            }

            continue


        # get the non-max frequent ones (the mutations)
        main_symbol = max(bps_freq, key=bps_freq.get)
        bps_freq.pop(main_symbol, None)
        bps_freq_norm = list(map(lambda x: x / bps, bps_freq.values()))
        bps_freq_norm = dict(zip(bps_freq.keys(), bps_freq_norm))

        mutation_locations = {}

        # for each symbol, get locations
        for symbol in bps_freq.keys():

            indices = [i for i in range(len(location_bps)) if location_bps[i] == symbol]
            mutation_locations[symbol] = indices


        mutations[bp] = {
            'bps_freq': bps_freq,
            'bps_freq_norm': bps_freq_norm,
            'mutation_locations': mutation_locations
        }

    return mutations

'''
# for each country, see how many mutations there are (only if there is a change)
for country in df_countries.keys():

    country_df = df_countries[country]

    # only look at the ones with 10 or more rows
    if country_df.shape[0] < 10:
        continue

    # calculate all the mutations
    m = mutations(country_df)

    # save one pickle per country
    with open(os.path.join('mutation_pickles', country), 'wb') as handle:
        pickle.dump(m, handle)
'''

df_dates = df.sort_values(by=['date'])
m = mutations(df_dates)
with open(os.path.join('mutation_pickles', 'all_dates'), 'wb') as handle:
    pickle.dump(m, handle)


# Look at the BDM main and density over time
# BDM rank order and color countries

'''
# list all the evolution files
evolution_folder = 'bdm_pickles/bdms_genes__data_covid_evolution'
files = os.listdir(evolution_folder)
files = list(filter(lambda x: not re.search('DS_Store', x), files))

# sequences are in this file
sequence_file = 'data_covid_evolution/gisaid_cov2020_sequences.fasta'

# read in the whole file
with open(sequence_file) as file:
    sequence_file = file.read()


# save all information to a dict
covid_evolution = {}


# get all the sequences loaded first

for file in files:

    data = pickle.load(open(os.path.join(evolution_folder, file), 'rb'))

    # find the sequence in the sequence file
    blob_head = re.sub('\|', '.{1}', file)
    blob = re.findall(blob_head + '[^\>]*', sequence_file)[0]

    # split by newline, join all except the first
    blob_parts = blob.split('\n')

    # get id and group, save to dict
    id = re.sub('>', '', blob_parts[0].split()[0])
    group = blob_parts[0].split()[-1]

    # make sequence from the other parts
    sequence = ''.join(blob_parts[1:])

    # throw away any small sequences
    if len(sequence) < 20000:
        continue

    try:
        date = re.findall('\d{4}\-\d{2}\-\d{2}', file)[0]
        date = datetime.strptime(date, '%Y-%m-%d')

        # throw away dates that are before 2019
        if date.year < 2019 or (date.year == 2019 and date.month < 11):
            continue

        country = file.split('|')[1]

    except:
        continue

    covid_evolution[file] = {
        'date': date,
        'country': country,
        'whole_bdm': data['whole_bdm'],
        'bdm_density': data['bdm_density']
    }


# make into a df to plot
df = pd.DataFrame.from_dict(covid_evolution)
df = df.transpose()


# plot results
plt.subplots(figsize=(10, 4))
plt.plot_date(x="date", y="whole_bdm", data=df, ms=2)
plt.xticks(rotation=70)
plt.xlabel('Date')
plt.ylabel('BDM')
plt.tight_layout()
plt.show()
'''


# Below is too hard because all the sequences are way too different from each other.
# For this script, we will load in each pre-computed bdm file for each "strain"
# load in each sequence, date, location
# Identify regions of change, mutations
# Get the BDM value for each region
# See how the BDM relates to region of change: Do regions with mutations have higher or lower BDM values?

'''
for file in files:

    data = pickle.load(open(os.path.join(evolution_folder, file), 'rb'))

    try:
        date = re.findall('\d{4}\-\d{2}\-\d{4}', file)[0]
        country = file.split('|')[1]
    except:
        continue

    # find the sequence in the sequence file
    blob_head = re.sub('\|', '.{1}', file)
    blob = re.findall(blob_head + '[^\>]*', sequence_file)[0]

    # split by newline, join all except the first
    blob_parts = blob.split('\n')

    # get id and group, save to dict
    id = re.sub('>', '', blob_parts[0].split()[0])
    group = blob_parts[0].split()[-1]

    # make sequence from the other parts
    sequence = ''.join(blob_parts[1:])

    # throw away any small sequences
    if len(sequence) < 20000:
        continue

    # N = Any of the ATCG
    # Trim away anything that isn't ATCG
    if re.search('[^ATCG]', sequence):
        sequence = re.sub('[^ATCG]', '', sequence)


    covid_evolution[file] = {
        'date': date,
        'country': country,
        'sequence': sequence
    }
'''
