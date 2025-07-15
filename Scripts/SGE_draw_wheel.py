import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import matplotlib.patches as mpatches
from matplotlib.patches import Arc
import math
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
    ax.set_aspect('equal')
    plt.show()
    #fig.show()
    return fig, ax


def main():
    seq = input("Enter the sequence: ")
    num_residues = len(seq)

    x_center, y_center, num_residues = generate_wheel_coordinates(num_residues)
    test = draw_wheel(seq, num_residues, x_center, y_center, legend=True)
    print(test)
    
main()
