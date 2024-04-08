import pandas as pd
import os
import math
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objs as go


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


fig = make_subplots()
for i, (processed_file, raw_file) in enumerate(matching_files):

    processed_file_path = os.path.join(processed_directory, processed_file)
    raw_file_path = os.path.join(raw_directory, raw_file)

    processed_df = pd.read_csv(processed_file_path)
    raw_df = pd.read_csv(raw_file_path)

    processed_df = processed_df[~processed_df['Reporter Identifier'].isin(excluded_probes) & processed_df['Reporter Identifier'].str.startswith('A_')]
    raw_df = raw_df[~raw_df['ProbeName'].isin(excluded_probes) & raw_df['ProbeName'].str.startswith('A_')]

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

   
    processed_df = processed_df.sort_values('VALUE', ascending=False)
    processed_df = processed_df.drop_duplicates(subset='Reporter Identifier', keep='first')
    top_20_reporters = processed_df.head(20)
    '''
    for _, row in top_20_reporters.iterrows():
        reporter = row['Reporter Identifier']
        value = row['VALUE']
        description = raw_df[raw_df['ProbeName'] == reporter]['Description'].iloc[0]
        print(f'Reporter: {reporter}, Description: {description}, Value: {value}')
    '''

    top_10_reporters_df = pd.DataFrame({
        'Reporter': top_20_reporters['Reporter Identifier'],
        'Description': [raw_df[raw_df['ProbeName'] == reporter]['Description'].iloc[0] for reporter in top_20_reporters['Reporter Identifier']],
        'Value': top_20_reporters['VALUE']
    })

    
    bar_chart = px.bar(top_10_reporters_df, x='Reporter', y='Value', hover_data=['Description']).data[0]
    fig.add_trace(bar_chart)

    
    button = dict(
        label = str(processed_file[:9]) + ' ' + factor_value_dict.get(processed_file, ''),
        method = 'update',
        args = [{'visible': [False]*len(matching_files)}]
    )

    button['args'][0]['visible'][i] = True
    buttons.append(button)

fig.update_layout(
    updatemenus=[
        go.layout.Updatemenu(
            buttons=buttons,
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.1,
            xanchor="left",
            y=1.1,
            yanchor="top"
        ),
    ],
    title='Top 20 Probes for Each Sample',
    xaxis_title='Probe',
    yaxis_title='Frequency Value'
)


fig.show()