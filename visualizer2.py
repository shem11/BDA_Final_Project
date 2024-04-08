import pandas as pd
import os
import math
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import plotly.io as py 
import numpy as np


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
 
    factor_value = factor_value_dict.get(file_name, "")
    print(factor_value)

 
    if factor_value not in ["2.5Gy whole-tissue", "0.5Gy 625-875 um", "2.5Gy 375-625 um", "2.5Gy 125-375 um"]:
        return "Ordinary"
    else:
        return "Radiated"

fig = make_subplots(rows=2, cols=1)

distance_mapping = {
    '2.5Gy whole-tissue': 0,  
    '0.5Gy 625-875 um': (625 + 875) / 2,
    '2.5Gy 375-625 um': (375 + 625) / 2,
    '2.5Gy 125-375 um': (125 + 375) / 2
}

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

    file_type = check_file_type(processed_file)

    merged_df['Type'] = file_type
    merged_df['Sample'] = str(processed_file)[:9] + ' ' + factor_value_dict.get(processed_file, "")
    merged_df['Distance'] = factor_value_dict.get(processed_file, "")

    merged_df = merged_df.drop_duplicates(subset='Reporter Identifier', keep='first')
    all_merged_dfs.append(merged_df)

all_data = pd.concat(all_merged_dfs)

ordinary_df = all_data[all_data['Type'] == 'Ordinary'].nlargest(200, 'VALUE')
radiated_df = all_data[all_data['Type'] != 'Ordinary'].nlargest(200, 'VALUE')
ordinary_df['Description_Sample'] = ordinary_df['Description'] + ' (Sample: ' + ordinary_df['Sample'] + ')'
radiated_df['Description_Sample'] = radiated_df['Description'] + ' (Sample: ' + radiated_df['Sample'] + ')'

distance_mapping = {
    '2.5Gy whole-tissue': 0,  
    '0.5Gy 625-875 um': (625 + 875) / 2,
    '2.5Gy 375-625 um': (375 + 625) / 2,
    '2.5Gy 125-375 um': (125 + 375) / 2
}

radiated_df['Distance'] = radiated_df['Sample'].map(distance_mapping)



ordinary_pivot = ordinary_df.pivot_table(index='Reporter Identifier', columns='Sample', values='VALUE')
radiated_pivot = radiated_df.pivot_table(index='Reporter Identifier', columns='Sample', values='VALUE')


ordinary_df['Type'] = 'Ordinary'
radiated_df['Type'] = 'Radiated'


df_combined = pd.concat([ordinary_df, radiated_df])
combined_pivot = df_combined.pivot_table(index='Reporter Identifier', columns='Sample', values='VALUE', aggfunc=np.sum)


hovertext_ordinary = ordinary_pivot.copy()
for i in range(len(hovertext_ordinary.index)):
    for j in range(len(hovertext_ordinary.columns)):
        hovertext_ordinary.iloc[i, j] = f'Sample: {hovertext_ordinary.columns[j]}<br>Probe: {hovertext_ordinary.index[i]}<br>Description: {ordinary_df.loc[ordinary_df["Reporter Identifier"] == hovertext_ordinary.index[i], "Description"].values[0]}<br>Value: {hovertext_ordinary.iloc[i, j]}'

hovertext_radiated = radiated_pivot.copy()
for i in range(len(hovertext_radiated.index)):
    for j in range(len(hovertext_radiated.columns)):
        hovertext_radiated.iloc[i, j] = f'Sample: {hovertext_radiated.columns[j]}<br>Probe: {hovertext_radiated.index[i]}<br>Description: {radiated_df.loc[radiated_df["Reporter Identifier"] == hovertext_radiated.index[i], "Description"].values[0]}<br>Value: {hovertext_radiated.iloc[i, j]}'

hovertext_combined = combined_pivot.copy()
for i in range(len(hovertext_combined.index)):
    for j in range(len(hovertext_combined.columns)):
        hovertext_combined.iloc[i, j] = f'Sample: {hovertext_combined.columns[j]}<br>Probe: {hovertext_combined.index[i]}<br>Description: {all_data.loc[all_data["Reporter Identifier"] == hovertext_combined.index[i], "Description"].values[0]}<br>Value: {hovertext_combined.iloc[i, j]}'



fig = go.Figure()
fig.add_trace(
    go.Heatmap(
        z=ordinary_pivot.values,
        x=ordinary_pivot.columns,
        y=ordinary_pivot.index,
        colorscale='Blues',
        name='Ordinary',
        showscale=True,
        colorbar=dict(
            titleside='right',
            tickmode='auto',
            nticks=10
        ),
        hovertext=hovertext_ordinary.values,
        hoverinfo='text',
        visible=True
    )
)

fig.add_trace(
    go.Heatmap(
        z=radiated_pivot.values,
        x=radiated_pivot.columns,
        y=radiated_pivot.index,
        colorscale='Reds',
        name='Radiated',
        showscale=True,
        colorbar=dict(
            titleside='right',
            tickmode='auto',
            nticks=10
        ),
        hovertext=hovertext_radiated.values,
        hoverinfo='text',
        visible=False 
    )
)

fig.add_trace(
    go.Heatmap(
        z=combined_pivot.values,
        x=combined_pivot.columns,
        y=combined_pivot.index,
        colorscale='Viridis',
        name='Combined',
        showscale=True,
        colorbar=dict(
            titleside='right',
            tickmode='auto',
            nticks=10
        ),
        hovertext=hovertext_combined.values,
        hoverinfo='text',  
        visible=False  
    )
)

fig.update_layout(
    updatemenus=[
        dict(
            buttons=list([
                dict(
                    args=[{"visible": [True, False, False]}, {"transition": {"duration": 500}}],
                    label="Ordinary",
                    method="restyle"
                ),
                dict(
                    args=[{"visible": [False, True, False]}, {"transition": {"duration": 500}}],
                    label="Radiated",
                    method="restyle"
                ),
                dict(
                    args=[{"visible": [False, False, True]}, {"transition": {"duration": 500}}],
                    label="Combined",
                    method="restyle"
                )
            ]),
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.1,
            xanchor="left",
            y=1.1,
            yanchor="top"
        ),
    ],
    xaxis_title="Tissue Sample",
    yaxis_title="Frequency Value",
    title="Heat-map of Frequency Values per Tissue Sample\n"
)


fig.show()