
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from pymol import cmd

bard1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20250508_BARD1scores_update_FILTERED.xlsx'
brca1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20240830_BRCA1_SGE_AllScores.xlsx'

bard1_data = pd.read_excel(bard1_file)
brca1_data = pd.read_excel(brca1_file)

#Helix residues for BARD1: [24, 48] and [98, 117]
#Helix residues for BRCA1: [7, 22] and [80, 97]
bard1_helix_residues = list(range(24, 49)) + list(range(98, 118))
brca1_helix_residues = list(range(7, 23)) + list(range(80, 98))

bard1_data = bard1_data.rename(columns = {'simplified_consequence': 'Consequence'})
brca1_data = brca1_data.rename(columns = {'snv_score_minmax': 'score'})
bard1_data = bard1_data.loc[bard1_data['Consequence'].isin(['missense_variant'])]
brca1_data = brca1_data.loc[brca1_data['Consequence'].isin(['missense_variant'])]

bard1_data['AApos'] = bard1_data['amino_acid_change'].transform(lambda x: x[1:-1])

brca1_data['AApos'] = brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][3:-3])

bard1_data = bard1_data.loc[bard1_data['AApos'].astype(int).isin(bard1_helix_residues)]
brca1_data = brca1_data.loc[brca1_data['AApos'].astype(int).isin(brca1_helix_residues)]


bard1_collapsed = bard1_data.groupby('AApos').agg({'score': 'median'}).reset_index()
brca1_collapsed = brca1_data.groupby('AApos').agg({'score': 'median'}).reset_index()

print(brca1_collapsed.head())