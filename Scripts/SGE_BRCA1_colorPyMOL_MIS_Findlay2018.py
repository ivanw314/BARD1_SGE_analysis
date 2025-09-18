from pymol import cmd
import pandas as pd
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt


#This script isn't super user-friendly but here's the run down:
#First, please make sure the correct file is selected
#Second, make sure you have the correct region selected, regions are dataset  
#specific and commented out to describe which domain and PDB ID it is from
#Third, 
#BE SURE TO SET THE OFFSET IN PYMOL and set coloring based on chains in the
#make_residue_values() function

file = '/Users/ivan/Desktop/20240830_BRCA1_SGE_AllScores.xlsx' #BRCA1 SGE scores file


#This block contains the list of tuples corresponding to which regions in the 
#data map to which structural domains in the provided PDB structure
#regions = [(1,301)] #BRCA1 RING from 1JM7
regions = [(4945,5565)] #BRCA1 BRCT from 1T29

def read_scores(file): #Reads and filters the score files
    excel = pd.read_excel(file) #Reads excel file into df
    data = excel[['target','pos','Consequence','snv_score_minmax']] #pulls out these relevant columns
    return data

def pull_scores(data, regions): #Filters for missense scores
    coords = [] #List to hold coordinates for each codon
    consequence = ['missense_variant'] #condition upon which data is fitlered (missense only)
    for start, end in regions: #loop that creates the coordinates for all codons
        for i in range(start, end + 1):
            coords.append(i)
    filtered = data[data['pos'].isin(coords) & data['Consequence'].isin(consequence)] #filters data for missense only and only within the specified coordinates
    
    filtered = filtered.reset_index(drop = True) #resets index to look nice

    codon_num = int(len(coords) / 3) #Gets number of codons
        
    return filtered, codon_num, coords

def make_residue_values(data, num, coords): #Gets mean score for all codons
    codon_score = {} #Dictionary that will hold scores for all codons in form of {Codon Number:Mean Score}
    codon_lists = [] #List to hold lists that contain the three position IDs for each codon
    
    #This loop creates the lists that are in the codon_lists element
    i = 0
    while (i) < len(coords): #for this while loop, have condition be i if starting region coordinate isn't 1, if it is 1, then have (i+3)
        codon = [] #empty list that will hold all the positions
        for i in range(i, i + 3): #for loop iterates through and creates the 3 positions
            codon.append(coords[i]) #appends to empty list
        codon_lists.append(codon) #appends to codon)list

        i += 1

    #This loop gets the mean score for each codon and makes a dictionary element    
    j = 0
    while j < num:
        data_filt = data.copy() #copy of data for each loop 
        data_filt = data[data['pos'].isin(codon_lists[j])] #data filtered so that you get data for just one codon
        mean = data_filt['snv_score_minmax'].min() #mean of scores is taken
        
        #offset needed to assign correct residue number in PyMOL structure
        offset = 1649 #offset in PyMOL structure (what is the number of the AA that starts the coloring?) (1 for RING, 1649 for BRCT)
        
        codon_score[j+ offset] = mean #mean score found and assigned to dictionary with key of base codon number + offset
        
        #this code can be used as a way to use medians for color mapping instead of mean
        #median = data_filt['snv_score_minmax'].median()
        #codon_score[j+1] = median
        j += 1
        
    return codon_score

def normalize_values(values): #Normalizes all values between 0 and 1 for coloring
    # First clamp all values between 0 and 1
    clamped_values = {k: min(max(v, -1.5), 0) for k, v in values.items()}
    clamped_values = {k: v for k, v in clamped_values.items() if not pd.isna(v)} #Filters out NA values
    
    # Get min and max of clamped values
    min_val = min(clamped_values.values())
    max_val = max(clamped_values.values())
    
    # Avoid division by zero if all values are the same
    if max_val == min_val:
        return {k: min_val for k in values.keys()}
        
    # Perform normalization
    return {k: (v - min_val) / (max_val - min_val) for k, v in clamped_values.items()}


colors = [(1, 1, 1), (1, 0, 0)]  # Gray -> Red
n_bins = 1000  # Number of bins for the colormap
cmap_name = 'gray_to_red'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins) #makes the color map 

def create_colorbar_legend():
    # Create figure with vertical proportions
    fig, ax = plt.subplots(figsize=(1.5, 5))
    fig.subplots_adjust(left=0.6)  # Adjust for labels on left
    
    # Reverse the colormap to match your get_color inversion
    reversed_cmap = custom_cmap.reversed()
    
    # Now use normal ordering since we reversed the colormap
    norm = plt.Normalize(vmin=-1.5, vmax=0)
    sm = cm.ScalarMappable(cmap=reversed_cmap, norm=norm)
    sm.set_array([])
    
    # Create vertical colorbar
    cbar = plt.colorbar(sm, cax=ax, orientation='vertical')
    
    # Move ticks and labels to the left
    cbar.ax.yaxis.set_ticks_position('left')
    cbar.ax.yaxis.set_label_position('left')
    
    # Labels now directly correspond to your data range
    cbar.set_ticks([-1.5, -1, -0.5, 0])
    cbar.set_label('Score', rotation=90, labelpad=20)  # Changed rotation to 90
    
    #plt.show()
    return fig

def get_color(value): #Gets color for each residue from mean score
    return custom_cmap(1 - value) #Inverts the colors, Gray is high score and Red is low


def main():
    data = read_scores(file) #Reads data
    filtered, num, coords = pull_scores(data, regions) #Gets filtered scores
    residue_values = make_residue_values(filtered, num, coords ) #Makes per-residue mean scores
    normalized_values = normalize_values(residue_values) #Scores normalized to between 0 and 1
    legend = create_colorbar_legend()
    #legend.savefig('/Users/ivan/Desktop/BARD1_draft_figs/fig5a_BRCA1_legend.png', dpi = 500)
    #this block does the coloring
    for residue, value in normalized_values.items(): 
        color_name = f'color_A_{residue}' #color_A specifies chain A
        color = get_color(value) #Gets color from color map
        cmd.set_color(color_name, [color[0], color[1], color[2]])  # RGB values
        cmd.color(color_name, f'chain A and resi {residue}') #chain  specifies chain A, change if not chain A
        #cmd.color(color_name, f'resi {residue}')
    cmd.show('cartoon')

main()