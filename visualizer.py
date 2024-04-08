import pandas as pd
import os
import math
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import plotly.io as py 


processed_directory = './Cleaned_processed_microarrays'
raw_directory = './Cleaned_raw_microarrays'
processed_filenames = os.listdir(processed_directory)
raw_filenames = os.listdir(raw_directory)
matching_files = [(processed_file, raw_file) for processed_file in processed_filenames for raw_file in raw_filenames if processed_file[:9] == raw_file[:9]]
excluded_probes = ['DarkCorner', '(-)3xSLv1', '(+)E1A_r60_a104', '(+)E1A_r60_a22', '(+)E1A_r60_a107', '(+)E1A_r60_a20', '(+)E1A_r60_n9', '(+)E1A_r60_n11', '(+)E1A_r60_3', '(+)E1A_r60_1']  # Some of the of the probes that are excluded are either controls or used for manipulation(rotation)

num_files = len(matching_files)
grid_size = math.ceil(math.sqrt(num_files))

gene_chip_df = pd.read_csv('GeneChip_assay_cleaned_new.csv')
factor_value_dict = dict(zip(gene_chip_df['Derived Array Data File'], gene_chip_df['Factor Value[Ionizing Radiation]']))

buttons = []
all_merged_dfs = []

def check_file_type(file_name):
    # Get the value from factor_value_dict
    factor_value = factor_value_dict.get(file_name, "")

    if factor_value not in ["2.5Gy whole-tissue", "0.5Gy 625-875 um", "2.5Gy 375-625 um", "2.5Gy 125-375 um"]:
        return "Ordinary"
    else:
        return "Radiated"

fig = make_subplots()

for i, (processed_file, raw_file) in enumerate(matching_files):
    processed_file_path = os.path.join(processed_directory, processed_file)
    raw_file_path = os.path.join(raw_directory, raw_file)


    processed_df = pd.read_csv(processed_file_path)
    raw_df = pd.read_csv(raw_file_path)
    processed_df = processed_df[~processed_df['Reporter Identifier'].isin(excluded_probes) & processed_df['Reporter Identifier'].str.startswith('A_')]
    raw_df = raw_df[~raw_df['ProbeName'].isin(excluded_probes) & raw_df['ProbeName'].str.startswith('A_')]

    processed_df = processed_df.drop_duplicates(subset='Reporter Identifier', keep='first')

    merged_df = pd.merge(processed_df, raw_df, left_on='Reporter Identifier', right_on='ProbeName')
    reporter_frequency = merged_df['Reporter Identifier'].value_counts()

    print(f'Sample: {str(processed_file)[:9]} {factor_value_dict.get(processed_file, "")}')
    '''
    # Print the 'Description' for the top 10 reporters
    top_10_reporters = reporter_frequency.nlargest(10)
    for reporter in top_10_reporters.index:
        description = raw_df[raw_df['ProbeName'] == reporter]['Description'].iloc[0]
        print(f'Reporter: {reporter}, Description: {description}')
    '''

    file_type = check_file_type(processed_file)

    merged_df['Type'] = file_type
    merged_df['Sample'] = str(processed_file)[:9]


    merged_df = merged_df.drop_duplicates(subset='Reporter Identifier', keep='first')
    all_merged_dfs.append(merged_df)

all_data = pd.concat(all_merged_dfs)

ordinary_df = all_data[all_data['Type'] == 'Ordinary'].nlargest(2000, 'VALUE')
radiated_df = all_data[all_data['Type'] != 'Ordinary'].nlargest(2000, 'VALUE')
ordinary_df['Description_Sample'] = ordinary_df['Description'] + ' (Sample: ' + ordinary_df['Sample'] + ')'
radiated_df['Description_Sample'] = radiated_df['Description'] + ' (Sample: ' + radiated_df['Sample'] + ')'

fig = go.Figure()
fig.add_trace(go.Scatter(
    x=ordinary_df['Reporter Identifier'], 
    y=ordinary_df['VALUE'], 
    mode='markers', 
    name='Ordinary',
    text=ordinary_df['Description_Sample']
))
fig.add_trace(go.Scatter(
    x=radiated_df['Reporter Identifier'], 
    y=radiated_df['VALUE'], 
    mode='markers', 
    name='Radiated',
    text=radiated_df['Description_Sample']
))

fig.update_layout(
    title='Top 2000 Probes',
    xaxis_title='Probe',
    yaxis_title='Frequency Value'
)

fig.show()