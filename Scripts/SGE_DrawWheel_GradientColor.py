import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

#Helical wheel code originally adapted from: Helixvis https://doi.org/10.21105/joss.01008 
# BARD1 and BRCA1 cutoffs for functional classification
bard1_cutoffs = [-3.438968 * 0.028675 + 0.009242,-3.018904 * 0.028675 + 0.009242]
brca1_cutoffs = [-1.328,-0.748]
# File paths
# Update these paths to the correct locations of your files
bard1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20250508_BARD1scores_update_FILTERED.xlsx'
brca1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20240830_BRCA1_SGE_AllScores.xlsx'

#Figure Saving Path
path = '/Users/ivan/Desktop/BARD1_draft_figs/4helix_helical_wheels/'
analysis_type = 'median_NP'  # Type of analysis. 'min', 'median', 'min_NP', 'median_NP' for minimum, mean score, or minimum/mean (proline substituions removed)
normalization_type = 'stdev' #Normalization type for color gradient. 'minmax' for min-max normalization based on 5th and 95th percentiles, 'stdev' for standard deviation normalization
pd.options.mode.chained_assignment = None

def read_process_data(bard1_file, brca1_file, type): #Reads and processes the data from the BARD1 and BRCA1 files

    aa_3to1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    } #Dictionary for converting 3-letter amino acid codes to 1-letter codes

    bard1_data = pd.read_excel(bard1_file) #Reads the BARD1 data from an Excel file
    brca1_data = pd.read_excel(brca1_file) #Reads the BRCA1 data from an Excel file

    #Helix residues for BARD1: [34, 48] and [98, 117]
    #Helix residues for BRCA1: [7, 22] and [80, 97]
    bard1_helix_1 = list(range(34, 49)) 
    bard1_helix_2 = list(range(99, 117))
    brca1_helix_1 = list(range(7, 23))
    brca1_helix_2 = list(range(80, 98))

    bard1_helix_residues = [bard1_helix_1, bard1_helix_2]
    brca1_helix_residues = [brca1_helix_1, brca1_helix_2]

    bard1_data = bard1_data.rename(columns = {'simplified_consequence': 'Consequence'}) #Renaming columns for consistency
    brca1_data = brca1_data.rename(columns = {'snv_score_minmax': 'score'}) #Renaming columns for consistency
    bard1_data = bard1_data.loc[bard1_data['Consequence'].isin(['missense_variant'])] #pulling only missense variants
    brca1_data = brca1_data.loc[brca1_data['Consequence'].isin(['missense_variant'])] #pulling only missense variants

    bard1_data['AApos'] = bard1_data['amino_acid_change'].transform(lambda x: x[1:-1]) #Gets amino acid position from the amino acid change
    bard1_data['amino_acid'] = bard1_data['amino_acid_change'].transform(lambda x: x[0]) #Gets original amino acid from the amino acid change
    bard1_data['amino_acid'] = bard1_data['amino_acid'] + bard1_data['AApos'] #Gets the full amino acid with position

    brca1_data['AApos'] = brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][3:-3]) #Gets amino acid position from the HGVS protein notation
    brca1_data['amino_acid'] = brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][0:3]) #Gets original amino acid from the HGVS protein notation
    brca1_data['amino_acid'] = brca1_data['amino_acid'].map(aa_3to1) #Remaps the 3-letter amino acid code to 1-letter code
    brca1_data['amino_acid'] = brca1_data['amino_acid'] + brca1_data['AApos'] #Gets the full amino acid with position

    bard1_mean = bard1_data['score'].mean() #Calculates the mean helix score for BARD1
    brca1_mean = brca1_data['score'].mean() #Calculates the mean helix score for BRCA1

    bard1_stdev =  bard1_data['score'].std()
    brca1_stdev = brca1_data['score'].std() #Calculates 2 standard deviations for BARD1 and BRCA1 helix scores

    bard1_min = bard1_data['score'].quantile(0.05)
    bard1_max = bard1_data['score'].quantile(0.95)

    brca1_min = brca1_data['score'].quantile(0.05)
    brca1_max = brca1_data['score'].quantile(0.95)
    final_dfs = {} #Dicitonary to store final dictionaries for each helix
    uncategorized_dfs = {} #Dictionary to store uncategorized dictionaries for each helix

    helix_list = ['helix_1', 'helix_2'] #List for iteration over helix names

    i = 0
    while i < len(bard1_helix_residues): #Iterates over the helix residues for BARD1 and collapses scores based on selected analysis
        elem = bard1_helix_residues[i] #Gets helix residues for the current helix
        data = bard1_data.loc[bard1_data['AApos'].astype(int).isin(elem)] #Gets data for the current helix residues

        #Aggregates the data based on the selected type of analysis
        if type == 'min':
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'median':
            to_return = data.groupby('AApos').agg({'score': 'median',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'min_NP':
            data['AAsub'] = data['amino_acid_change'].transform(lambda x: x[-1])
            data = data.loc[data['AAsub'] != 'P']
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'median_NP':
            data['AAsub'] = data['amino_acid_change'].transform(lambda x: x[-1])
            data = data.loc[data['AAsub'] != 'P']
            to_return = data.groupby('AApos').agg({'score': 'median',
                                                   'amino_acid': 'first'}).reset_index()
        else:
            raise ValueError("Invalid type specified. Choose from 'min', 'median', 'min_NP', 'median_NP'.")

        to_return['median_consequence'] = 1 #Sets default function to 1 (indeterminate) (color to yellow)
        to_return.loc[to_return['score'] <= bard1_cutoffs[0], 'median_consequence'] = 3 #Sets low score variants to abnormal color to red
        to_return.loc[to_return['score'] >= bard1_cutoffs[1], 'median_consequence'] = 2 #Sets high score variatns to normal, color to blue
        to_return['AApos'] = to_return['AApos'].astype(int) #Converts amino acid position to integer for sorting
        to_return.sort_values(by='AApos', inplace=True) #Sorts in ascending order by amino acid position

        if helix_list[i] ==  'helix_2': #2nd helix should be sorted in descending order
            to_return.sort_values(by='AApos', ascending=False, inplace=True)

        to_return_dict = dict(zip(to_return['amino_acid'], to_return['median_consequence'])) #Zips and returns a dictionary with amino acid as key and median consequence as value

        final_dfs['bard1_' + helix_list[i]] = to_return_dict #Adds the dictionary to the final_dfs dictionary with key as 'bard1_' + helix name
        uncategorized_dfs['bard1_' + helix_list[i]] = to_return #Adds the uncategorized dictionary to the uncategorized_dfs dictionary with key as 'bard1_' + helix name
        i += 1
    
    i = 0 
    while i < len(brca1_helix_residues): #Analogous process for BRCA1 helix residues
        elem = brca1_helix_residues[i]
        data = brca1_data.loc[brca1_data['AApos'].astype(int).isin(elem)]

        if type == 'min':
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'median':
            to_return = data.groupby('AApos').agg({'score': 'median',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'min_NP':
            data['AAsub'] = data['hgvs_pro'].transform(lambda x: x[-3:])
            data = data.loc[data['AAsub'] != 'Pro']
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'median_NP':
            data['AAsub'] = data['hgvs_pro'].transform(lambda x: x[-3:])
            data = data.loc[data['AAsub'] != 'Pro']
            to_return = data.groupby('AApos').agg({'score': 'median',
                                                   'amino_acid': 'first'}).reset_index()
        else:
            raise ValueError("Invalid type specified. Choose from 'min', 'median', 'min_NP', 'median_NP'.")

        to_return['median_consequence'] = 1
        to_return.loc[to_return['score'] <= brca1_cutoffs[0], 'median_consequence'] = 3
        to_return.loc[to_return['score'] >= brca1_cutoffs[1], 'median_consequence'] = 2
        to_return['AApos'] = to_return['AApos'].astype(int)

        to_return.sort_values(by='AApos', inplace=True)
        if helix_list[i] ==  'helix_2':
            to_return.sort_values(by='AApos', ascending=False, inplace=True)
        to_return_dict = dict(zip(to_return['amino_acid'], to_return['median_consequence']))

        stats_dict = {'bard1': (bard1_mean, bard1_stdev, bard1_min, bard1_max), 'brca1': (brca1_mean, brca1_stdev, brca1_min, brca1_max)} #Dictionary to hold the statistics for each gene
        final_dfs['brca1_' + helix_list[i]] = to_return_dict    
        uncategorized_dfs['brca1_' + helix_list[i]] = to_return
        i += 1

    dict_keys = list(final_dfs.keys())

    return final_dfs,uncategorized_dfs,stats_dict, dict_keys

def gradient_color(helical_dict, stats_dict, norm_type):
    colors = [(1, 1, 1.0), (1, 0, 0)]  # White -> Red
    n_bins = 1000  # Number of bins for the colormap
    cmap_name = 'gray_to_red'
    custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins) #makes the color map 

    #custom_cmap = plt.cm.Reds
    stats_keys = list(stats_dict.keys())
    helix_df_keys = list(helical_dict.keys())
    to_color = {}
    for gene in stats_keys: #Iterates over the genes (BARD1 and BRCA1)
        mean, stdev, min, max = stats_dict[gene]

        if norm_type == 'minmax':
            lwrbound = min
            uprbound = max
        elif norm_type == 'stdev':
            lwrbound = mean - 1*stdev #Calculates lower bound for the gradient color
            uprbound = mean + 1*stdev #Calculates upper bound for the gradient color

        for helix in helix_df_keys: #Iterates over the helix names
            if gene in helix:
                helix_df = helical_dict[helix] #Gets the helix dictionary for the current helix
                if norm_type == 'stdev':
                    helix_df['normalized_score'] = np.interp(helix_df['score'], (lwrbound, uprbound), (0, 1)) #Normalizes the score for the gradient color
                elif norm_type == 'minmax':
                    print(lwrbound, uprbound)
                    helix_df['normalized_score'] = helix_df['score'].apply(lambda x: (x - lwrbound) / (uprbound - lwrbound)) #Normalizes the score for the gradient color
                print(helix_df)
                helix_df['color'] = helix_df['normalized_score'].apply(lambda x: custom_cmap(1-x))
                to_color[helix] = dict(zip(helix_df['amino_acid'], helix_df['color'])) #Creates a dictionary with amino acid as key and normalized score as value
            else:
                continue

    #print(to_color)

    return to_color
def generate_wheel_coordinates(dict): #Generates coordinates for the helical wheel plot based on length of sequence

    n = len(dict)  #Number of residues from length of the dictionary 
    x_center = []
    y_center = []
    num_residues = n  # Maximum number of residues
    for i in range(n):
        # Start at 90° (top) and go clockwise (-100° each step)
        angle = (90 - i * 100) * np.pi / 180
        x = np.cos(angle) * 0.85  # 0.85 is the radius from the original
        y = np.sin(angle) * 0.85
        x_center.append(x)
        y_center.append(y)
    return np.array(x_center), np.array(y_center), num_residues

def missense_draw_wheel(sequence, path, resi_dict, helix_name, num_residues, x_array,  y_array, colors = ["gray", "yellow", "blue", "red"], labels = False, labelcolor = "black", legend = False): #Draws the helical wheel plot for a given sequence
    "draw helix"
    min_num = 2 # Minimum number of residues
    max_num = num_residues # Maximum number of residues
    num_colors = 4
    num_resid = len(sequence) # Length of the sequence

    residues = resi_dict # Dictionary of residues with their corresponding colors based on aggregated functional consequence
    residue_keys = list(residues.keys()) # List of residue keys

    colors = []
    for key in residue_keys:
        colors.append(residues[key]) # Appends the colors for each residue based on the functional consequence

    '''
    #Error checking
    if num_resid not in range(min_num, max_num + 1):
        return "ERROR: sequence must have between 2 and 18 (inclusive) characters."
    if len(colors) != 4:
        return "ERROR: parameter `colors` has missing or too many colors."
    for i in range(len(colors)):
        if colors[i] not in pltcol.cnames:
            return "ERROR: parameter `colors` has invalid colors." 

    '''

    #Builds helical wheel
    x_center = x_array
    y_center = y_array
    x_center = x_center/2 + 0.5
    y_center = y_center/2 + 0.5
    circle_radius = 0.0725
    circle_data = pd.DataFrame(data={'x': x_center[0:num_resid], 
        'y': y_center[0:num_resid], 'type': range(num_resid)}
        )
    circle_data['color'] = pd.Series(dtype = 'object')
    for i in range(num_resid):
        if sequence[i] not in residues:
            return "ERROR: " + sequence[i] + " is not a valid one-letter code for an amino acid."
        
        circle_data['type'][i] = 0
        
        circle_data.at[i, 'color'] = colors[i]
        #circle_data['color'][i] = 3
        
    segment_data = pd.DataFrame(data={'xstart': x_center[0:num_resid - 1], 
        'ystart': y_center[0:num_resid - 1], 'xend': x_center[1:num_resid], 
        'yend': y_center[1:num_resid]})
    
    fig, ax = plt.subplots()
    for i in range(num_resid - 1):
        plt.plot([segment_data['xstart'][i], segment_data['xend'][i]], [segment_data['ystart'][i], segment_data['yend'][i]], 'ro-', color = 'black')
        
    for i in range(num_resid):
        circle = plt.Circle((circle_data['x'][i], circle_data['y'][i]), circle_radius, clip_on = False, zorder = 10, facecolor= circle_data['color'][i], edgecolor = 'black')
        ax.add_artist(circle)
        if labels:
            if helix_name == 'bard1_helix_1':
                ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=10, rotation = 120, ha="center", va = "center", color = labelcolor)
            elif helix_name == 'bard1_helix_2' or helix_name == 'brca1_helix_1':
                ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=10, rotation = -90, ha="center", va = "center", color = labelcolor)
            elif helix_name == 'brca1_helix_2':
                ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=10, rotation = 180, ha="center", va = "center", color = labelcolor)
    if legend:
        restypes = set(circle_data['type'])
        handleid = []
        nonpolar = mpatches.Patch(color = colors[0], label = 'hydrophobic')
        polar = mpatches.Patch(color = colors[1], label = 'polar')
        basic = mpatches.Patch(color = colors[2], label = 'basic')
        acidic = mpatches.Patch(color = colors[3], label = 'acidic')
        if 0 in restypes:
            handleid = [nonpolar]
            
        if 1 in restypes:
            if bool(handleid):
                handleid.append(polar)
            else:
                handleid = [polar]
                
        if 2 in restypes:
            if bool(handleid):
                handleid.append(basic)
            else:
                handleid = [basic]
        
        if 3 in restypes:
            if bool(handleid):
                handleid.append(acidic)
            else:
                handleid = [acidic]
                
        plt.legend(handles = handleid, loc='center left', bbox_to_anchor=(1.04, 0.5))
        
    plt.axis('off')
    plt.title(helix_name, 
              fontsize=12,
              y = 1.05)
    ax.set_aspect('equal')
    plt.show()
    #fig.savefig(path + helix_name + '.png', bbox_inches='tight', dpi=500, transparent=True)
    #fig.show()
    return fig, ax

def main():
    helical_dicts, uncategorized_dicts, stats, helices = read_process_data(bard1_file, brca1_file, type = analysis_type)
    to_color = gradient_color(uncategorized_dicts, stats, norm_type = normalization_type)


    for helix in helices:
        helix_dict = to_color[helix]
        x_center, y_center, num_residues = generate_wheel_coordinates(helix_dict)
        seq_list = list(helix_dict.keys())
        #print(helix, helix_dict)
        missense_draw_wheel(seq_list, path, helix_dict, helix, num_residues, x_center, y_center, labels = True)


main()
