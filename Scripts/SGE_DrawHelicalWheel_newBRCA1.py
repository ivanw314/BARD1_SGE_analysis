import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import matplotlib.patches as mpatches

#Helical wheel code originally adapted from: Helixvis https://doi.org/10.21105/joss.01008 
# BARD1 and BRCA1 cutoffs for functional classification

brca1_cutoffs = [-1.328,-0.748]
# File paths
# Update these paths to the correct locations of your files
bard1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20250825_BARD1snvscores_filtered.xlsx'
brca1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20240830_BRCA1_SGE_AllScores.xlsx'

#Figure Saving Path
path = '/Users/ivan/Desktop/BARD1_draft_figs/'
analysis_type = 'min_NP'  # Type of analysis. 'min', 'mean', 'min_NP', 'mean_NP' for minimum, mean score, or minimum/mean (proline substituions removed)
pd.options.mode.chained_assignment = None

def get_bard1_thresholds(bard1_file):
        # find the GMM thresholds
    target_value = 0.950
    bard1df = pd.read_excel(bard1_file)
    # Calculate the absolute difference for the Normal (N) density
    diffN = (bard1df['gmm_density_normal'] - target_value).abs()
    # Find the index of the minimum difference
    closest_index = diffN.idxmin()
    # Retrieve the row with the closest value
    closest_row_n = bard1df.loc[closest_index]

    # now repeat that for the abnormal density
    # Calculate the absolute difference
    diffA = (bard1df['gmm_density_abnormal'] - target_value).abs()
    # Find the index of the minimum difference
    closest_index = diffA.idxmin()
    # Retrieve the row with the closest value
    closest_row_a = bard1df.loc[closest_index]

    # now we get the scores that are the closest to the (n)ormal and (a)bnormal thresholds
    score_n_95 = closest_row_n['score']
    score_a_95 = closest_row_a['score']

    bard1_cutoffs = [score_n_95, score_a_95]
    return bard1_cutoffs


def read_process_data(bard1_file, brca1_file, bard1_cutoffs, type): #Reads and processes the data from the BARD1 and BRCA1 files

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

    bard1_data = bard1_data.rename(columns = {'consequence': 'Consequence'}) #Renaming columns for consistency
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

    final_dfs = {} #Dicitonary to store final dictionaries for each helix
    helix_list = ['helix_1', 'helix_2'] #List for iteration over helix names

    i = 0
    while i < len(bard1_helix_residues): #Iterates over the helix residues for BARD1 and collapses scores based on selected analysis
        elem = bard1_helix_residues[i] #Gets helix residues for the current helix
        data = bard1_data.loc[bard1_data['AApos'].astype(int).isin(elem)] #Gets data for the current helix residues

        #Aggregates the data based on the selected type of analysis
        if type == 'min':
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'mean':
            to_return = data.groupby('AApos').agg({'score': 'mean',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'min_NP':
            data['AAsub'] = data['amino_acid_change'].transform(lambda x: x[-1])
            data = data.loc[data['AAsub'] != 'P']
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'mean_NP':
            data['AAsub'] = data['amino_acid_change'].transform(lambda x: x[-1])
            data = data.loc[data['AAsub'] != 'P']
            to_return = data.groupby('AApos').agg({'score': 'mean',
                                                   'amino_acid': 'first'}).reset_index()

        to_return['median_consequence'] = 1 #Sets default function to 1 (indeterminate) (color to yellow)
        to_return.loc[to_return['score'] <= bard1_cutoffs[0], 'median_consequence'] = 3 #Sets low score variants to abnormal color to red
        to_return.loc[to_return['score'] >= bard1_cutoffs[1], 'median_consequence'] = 2 #Sets high score variatns to normal, color to blue
        to_return['AApos'] = to_return['AApos'].astype(int) #Converts amino acid position to integer for sorting
        to_return.sort_values(by='AApos', inplace=True) #Sorts in ascending order by amino acid position

        if helix_list[i] ==  'helix_2': #2nd helix should be sorted in descending order
            to_return.sort_values(by='AApos', ascending=False, inplace=True)

        to_return = dict(zip(to_return['amino_acid'], to_return['median_consequence'])) #Zips and returns a dictionary with amino acid as key and median consequence as value
        final_dfs['bard1_' + helix_list[i]] = to_return #Adds the dictionary to the final_dfs dictionary with key as 'bard1_' + helix name
        i += 1
    
    i = 0 
    while i < len(brca1_helix_residues): #Analogous process for BRCA1 helix residues
        elem = brca1_helix_residues[i]
        data = brca1_data.loc[brca1_data['AApos'].astype(int).isin(elem)]

        if type == 'min':
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'mean':
            to_return = data.groupby('AApos').agg({'score': 'mean',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'min_NP':
            data['AAsub'] = data['hgvs_pro'].transform(lambda x: x[-3:])
            data = data.loc[data['AAsub'] != 'Pro']
            to_return = data.groupby('AApos').agg({'score': 'min',
                                                   'amino_acid': 'first'}).reset_index()
        elif type == 'mean_NP':
            data['AAsub'] = data['hgvs_pro'].transform(lambda x: x[-3:])
            data = data.loc[data['AAsub'] != 'Pro']
            to_return = data.groupby('AApos').agg({'score': 'mean',
                                                   'amino_acid': 'first'}).reset_index()

        to_return['median_consequence'] = 1
        to_return.loc[to_return['score'] <= brca1_cutoffs[0], 'median_consequence'] = 3
        to_return.loc[to_return['score'] >= brca1_cutoffs[1], 'median_consequence'] = 2
        to_return['AApos'] = to_return['AApos'].astype(int)

        to_return.sort_values(by='AApos', inplace=True)
        if helix_list[i] ==  'helix_2':
            to_return.sort_values(by='AApos', ascending=False, inplace=True)
        to_return = dict(zip(to_return['amino_acid'], to_return['median_consequence']))

        final_dfs['brca1_' + helix_list[i]] = to_return
        i += 1

    dict_keys = list(final_dfs.keys())
    return final_dfs, dict_keys

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

def missense_draw_wheel(sequence, path, resi_dict, helix_name, num_residues, x_array,  y_array, colors = ["gray", "yellow", "white", "red"], labels = False, labelcolor = "black", legend = False): #Draws the helical wheel plot for a given sequence
    "draw helix"
    min_num = 2 # Minimum number of residues
    max_num = num_residues # Maximum number of residues
    num_colors = 4
    num_resid = len(sequence) # Length of the sequence

    residues = resi_dict # Dictionary of residues with their corresponding colors based on aggregated functional consequence
        
    #Error checking
    if num_resid not in range(min_num, max_num + 1):
        return "ERROR: sequence must have between 2 and 18 (inclusive) characters."
    if len(colors) != 4:
        return "ERROR: parameter `colors` has missing or too many colors."
    for i in range(len(colors)):
        if colors[i] not in pltcol.cnames:
            return "ERROR: parameter `colors` has invalid colors." 

    #Builds helical wheel
    x_center = x_array
    y_center = y_array
    x_center = x_center/2 + 0.5
    y_center = y_center/2 + 0.5
    circle_radius = 0.0725
    circle_data = pd.DataFrame(data={'x': x_center[0:num_resid], 
        'y': y_center[0:num_resid], 'color': range(num_resid), 'type': range(num_resid)})

    for i in range(num_resid):
        if sequence[i] not in residues:
            return "ERROR: " + sequence[i] + " is not a valid one-letter code for an amino acid."
        circle_data['color'][i] = colors[residues[sequence[i]]]
        circle_data['type'][i] = residues[sequence[i]]
  
    segment_data = pd.DataFrame(data={'xstart': x_center[0:num_resid - 1], 
        'ystart': y_center[0:num_resid - 1], 'xend': x_center[1:num_resid], 
        'yend': y_center[1:num_resid]})
    
    fig, ax = plt.subplots()
    for i in range(num_resid - 1):
        plt.plot([segment_data['xstart'][i], segment_data['xend'][i]], [segment_data['ystart'][i], segment_data['yend'][i]], 'ro-', color = 'black')
        
    for i in range(num_resid):
        circle = plt.Circle((circle_data['x'][i], circle_data['y'][i]), circle_radius, clip_on = False, zorder = 10, facecolor=circle_data['color'][i], edgecolor = 'black')
        ax.add_artist(circle)
        if labels:
            if helix_name == 'bard1_helix_1':
                ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=12, weight = 'bold', rotation = 120, ha="center", va = "center", color = labelcolor)
            elif helix_name == 'bard1_helix_2' or helix_name == 'brca1_helix_1':
                ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=12, weight = 'bold', rotation = -90, ha="center", va = "center", color = labelcolor)
            elif helix_name == 'brca1_helix_2':
                ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=12, weight = 'bold', rotation = 180, ha="center", va = "center", color = labelcolor)
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
    fig.savefig(path + 'fig_5b_' + helix_name + '.png', bbox_inches='tight', dpi=500, transparent=True)
    #fig.show()
    return fig, ax

def main():
    bard1_cutoffs = get_bard1_thresholds(bard1_file)
    helical_dicts, helices = read_process_data(bard1_file, brca1_file, bard1_cutoffs, type = analysis_type)
    
    for helix in helices:
        helix_dict = helical_dicts[helix]
        x_center, y_center, num_residues = generate_wheel_coordinates(helical_dicts[helix])
        seq_list = list(helical_dicts[helix].keys())
        print(helix, helix_dict)
        missense_draw_wheel(seq_list, path, helix_dict, helix, num_residues, x_center, y_center, labels = True)


main()
