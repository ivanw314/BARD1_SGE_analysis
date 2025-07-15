import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import matplotlib.patches as mpatches
from matplotlib.patches import Arc
import math

bard1_cutoffs = [-3.438968 * 0.028675 + 0.009242,-3.018904 * 0.028675 + 0.009242]
brca1_cutoffs = [-1.328,-0.748]
bard1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20250508_BARD1scores_update_FILTERED.xlsx'
brca1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20240830_BRCA1_SGE_AllScores.xlsx'
type = 'min_NP'  # 'min', 'mean', 'min_NP', 'mean_NP' for minimum, mean score, or minimum/mean (proline substituions removed)
pd.options.mode.chained_assignment = None

def generate_wheel_coordinates(n):
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

def draw_wheel(sequence, num_residues, x_array,  y_array, colors = ["gray", "yellow", "blue", "red"], labels = False, labelcolor = "black", legend = False):
    "draw helix"
    min_num = 2
    max_num = num_residues
    num_colors = 4
    num_resid = len(sequence)

    # 0 = hydrophobic, 1 = polar, 2 = basic, 3 = acidic
    residues = {"A":0, "R":2, "N":1, "D":3, "C":1,
                  "Q":1, "E":3, "G":0, "H":2, "I":0,
                  "L":0, "K":2, "M":0, "F":0, "P":0,
                  "S":1, "T":1, "W":0, "Y":0, "V":0}
        
    if num_resid not in range(min_num, max_num + 1):
        return "ERROR: sequence must have between 2 and 18 (inclusive) characters."
    if len(colors) != 4:
        return "ERROR: parameter `colors` has missing or too many colors."
    for i in range(len(colors)):
        if colors[i] not in pltcol.cnames:
            return "ERROR: parameter `colors` has invalid colors." 

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
            ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=10, ha="center", va = "center", color = labelcolor)
    print(circle_data)
        
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
    #plt.title(''''Amino Acid Wheel Representation'''', fontsize=12)
    ax.set_aspect('equal')
    plt.show()
    #fig.show()
    return fig, ax


def read_process_data(bard1_file, brca1_file, type):

    aa_3to1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
    }

    bard1_data = pd.read_excel(bard1_file)
    brca1_data = pd.read_excel(brca1_file)

    #Helix residues for BARD1: [34, 48] and [98, 117]
    #Helix residues for BRCA1: [7, 22] and [80, 97]
    bard1_helix_1 = list(range(34, 49)) 
    bard1_helix_2 = list(range(99, 117))
    brca1_helix_1 = list(range(7, 23))
    brca1_helix_2 = list(range(80, 98))

    bard1_helix_residues = [bard1_helix_1, bard1_helix_2]
    brca1_helix_residues = [brca1_helix_1, brca1_helix_2]

    bard1_data = bard1_data.rename(columns = {'simplified_consequence': 'Consequence'})
    brca1_data = brca1_data.rename(columns = {'snv_score_minmax': 'score'})
    bard1_data = bard1_data.loc[bard1_data['Consequence'].isin(['missense_variant'])]
    brca1_data = brca1_data.loc[brca1_data['Consequence'].isin(['missense_variant'])]

    bard1_data['AApos'] = bard1_data['amino_acid_change'].transform(lambda x: x[1:-1])
    bard1_data['amino_acid'] = bard1_data['amino_acid_change'].transform(lambda x: x[0])
    bard1_data['amino_acid'] = bard1_data['amino_acid'] + bard1_data['AApos']

    brca1_data['AApos'] = brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][3:-3])
    brca1_data['amino_acid'] = brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][0:3])
    brca1_data['amino_acid'] = brca1_data['amino_acid'].map(aa_3to1)
    brca1_data['amino_acid'] = brca1_data['amino_acid'] + brca1_data['AApos']

    final_dfs = {}
    helix_list = ['helix_1', 'helix_2']

    i = 0
    while i < len(bard1_helix_residues):
        elem = bard1_helix_residues[i]
        data = bard1_data.loc[bard1_data['AApos'].astype(int).isin(elem)]

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

        to_return['median_consequence'] = 1
        to_return.loc[to_return['score'] <= bard1_cutoffs[0], 'median_consequence'] = 3
        to_return.loc[to_return['score'] >= bard1_cutoffs[1], 'median_consequence'] = 2
        to_return['AApos'] = to_return['AApos'].astype(int)
        to_return.sort_values(by='AApos', inplace=True)

        if helix_list[i] ==  'helix_2':
            to_return.sort_values(by='AApos', ascending=False, inplace=True)

        to_return = dict(zip(to_return['amino_acid'], to_return['median_consequence']))
        final_dfs['bard1_' + helix_list[i]] = to_return
        i += 1
    
    i = 0 
    while i < len(brca1_helix_residues):
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


    print(final_dfs)

    return final_dfs

def missense_draw_wheel(sequence, num_residues, x_array,  y_array, colors = ["gray", "yellow", "blue", "red"], labels = False, labelcolor = "black", legend = False):
    "draw helix"
    min_num = 2
    max_num = num_residues
    num_colors = 4
    num_resid = len(sequence)

    # 0 = hydrophobic, 1 = polar, 2 = basic, 3 = acidic
    residues = {"A":0, "R":2, "N":1, "D":3, "C":1,
                  "Q":1, "E":3, "G":0, "H":2, "I":0,
                  "L":0, "K":2, "M":0, "F":0, "P":0,
                  "S":1, "T":1, "W":0, "Y":0, "V":0}
        
    if num_resid not in range(min_num, max_num + 1):
        return "ERROR: sequence must have between 2 and 18 (inclusive) characters."
    if len(colors) != 4:
        return "ERROR: parameter `colors` has missing or too many colors."
    for i in range(len(colors)):
        if colors[i] not in pltcol.cnames:
            return "ERROR: parameter `colors` has invalid colors." 

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
            ax.annotate(sequence[i], xy=(circle_data['x'][i], circle_data['y'][i]), zorder = 15, fontsize=10, ha="center", va = "center", color = labelcolor)
    print(circle_data)
        
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
    #plt.title(''''Amino Acid Wheel Representation'''', fontsize=12)
    ax.set_aspect('equal')
    plt.show()
    #fig.show()
    return fig, ax

def main():
    helical_dicts = read_process_data(bard1_file, brca1_file, type = 'min_NP')
    seq = input("Enter the sequence: ")
    num_residues = len(seq)

    x_center, y_center, num_residues = generate_wheel_coordinates(num_residues)
    error_return = draw_wheel(seq, num_residues, x_center, y_center, legend=False, labels=True)
    print(error_return)

main()
