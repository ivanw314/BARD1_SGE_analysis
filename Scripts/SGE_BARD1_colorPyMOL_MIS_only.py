from pymol import cmd
import pandas as pd
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

#This script generates colored ribbon structures in PyMOL based on SGE scores. Used in figure 6A. Can be used to build similar structures for all other domains

#This script isn't super user-friendly but here's the run down:
#Simply enter the region name and analysis type in the User-provided inputs block
#Make sure you have the correct datafile selected
#and set coloring based on chains in the make_residue_values() function


#User-provided inputs
region = 'RING' #Hardcode the region name here (RING, ARD, BRCT)
chain = 'B' #Chain in the PDB structure to color (A, B, etc)
analysis = 'min' #mininum or mean score used for coloring (min, mean)
file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/final_tables/BARD1_SGE_final_table.xlsx' #SGE Score file

show_legend = False #Whether to show the legend figure
save_legend = True #Whether to save the legend figure
#This block contains the list of tuples corresponding to which regions in the 
#data map to which structural domains in the provided PDB structure
ring = [(214809412,214809494),(214797061, 214797117), (214792298,214792445)] #RING (1JM7)
ard = [(214780560,214780601),(214769232,214769312),(214767482,214767654),(214752486,214752555)] # ARD (3C5R) 
brct = [(214745722,214745830),(214745067,214745159),(214730411,214730508),(214728685,214729008)] # BRCT (3FA2)

#Region Offsets (What amino acid residue does the structure start in PyMOL?)
region_offsets = {
    'RING': 26,
    'ARD': 425,
    'BRCT': 568
}

region_end = {'RING': 122, 'ARD': 545, 'BRCT': 777} #Last residue in each domain

expected_residues = list(range(region_offsets[region], region_end[region] + 1))

def read_scores(file): #Reads and filters the score files
    excel = pd.read_excel(file, sheet_name = 'scores') #Reads excel file into df
    excel = excel.loc[excel['var_type'].isin(['snv'])] #Filters for SNVs only
    excel = excel.loc[~excel['variant_qc_flag'].isin(['WARN'])]
    excel = excel.rename(columns = {'consequence': 'Consequence', 'score': 'snv_score'})
    data = excel[['exon','pos','Consequence','snv_score']] #pulls out these relevant columns
    return data

def get_region_info(region): #Gets the respective coordinates and offset for each domain
    
    if region == 'RING':
        region_coords = ring
        offset = region_offsets['RING']
        
    elif region == 'ARD':
        region_coords = ard
        offset = region_offsets['ARD']
        
    elif region == 'BRCT':
        region_coords = brct
        offset = region_offsets['BRCT']

    return region_coords, offset

def pull_scores(data, regions): #Filters for missense scores
    coords = [] #List to hold coordinates for each codon
    for start, end in regions: #loop that creates the coordinates for all codons
        for i in range(start, end + 1):
            coords.append(i)
    coords.sort(reverse = True)
    #filters data for missense only and only within the specified coordinates. Will include multi-consequence variants
    filtered = data[data['pos'].isin(coords) & data['Consequence'].str.contains('missense')] 
    #print(filtered)
    filtered = filtered.reset_index(drop = True) #resets index to look nice

    codon_num = int(len(coords) / 3) #Gets number of codons
        
    return filtered, codon_num, coords

def make_residue_values(data, num, coords, region_offset, analysis):
    codon_score = {}
    codon_lists = []
    
    #This loop creates the lists that are in the codon_lists element
    i = 0
    while (i) < len(coords): #for this while loop, have condition be i if starting region coordinate isn't 1, if it is 1, then have (i+3)
        codon = [] #empty list that will hold all the positions
        for i in range(i, i + 3): #for loop iterates through and creates the 3 positions
            codon.append(coords[i]) #appends to empty list
        codon_lists.append(codon) #appends to codon)list

        i += 1

    #This loop gets the min/mean score for each codon and makes a dictionary element    
    j = 0
    while j < num:
        data_filt = data.copy() #copy of data for each loop 
        data_filt = data[data['pos'].isin(codon_lists[j])] #data filtered so that you get data for just one codon
        
        if analysis == 'min':
            min_score = data_filt['snv_score'].min() #min of scores is taken
            
            #offset needed to assign correct residue number in PyMOL structure
            offset = region_offset #offset in PyMOL structure 
            
            codon_score[j+ offset] = min_score #min score found and assigned to dictionary with key of base codon number + offset
        
        elif analysis == 'mean':
            mean_score = data_filt['snv_score'].mean() #mean of scores is taken
            
            #offset needed to assign correct residue number in PyMOL structure
            offset = region_offset #offset in PyMOL structure 
            
            codon_score[j+ offset] = mean_score #mean score found and assigned to dictionary with key of base codon number + offset
 
        j += 1
        
    return codon_score

def normalize_values(values): #Normalizes all values between 0 and 1 for coloring
    # First clamp all values between 0 and 1
    clamped_values = {k: min(max(v, -0.2), 0) for k, v in values.items()}
    clamped_values = {k: v for k, v in clamped_values.items() if not pd.isna(v)} #Filters out NA values
    
    # Get min and max of clamped values
    min_val = min(clamped_values.values())
    max_val = max(clamped_values.values())
    
    # Avoid division by zero if all values are the same
    if max_val == min_val:
        return {k: min_val for k in values.keys()}
        
    # Perform normalization
    return {k: (v - min_val) / (max_val - min_val) for k, v in clamped_values.items()}


colors = [(1.0, 1.0, 1.0), (1, 0, 0)]  # White -> Red
n_bins = 1000  # Number of bins for the colormap
cmap_name = 'gray_to_red'
custom_cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins) #makes the color map 

def create_colorbar_legend():
    # Create figure with horizontal proportions
    fig, ax = plt.subplots(figsize=(1, 0.5))
    fig.subplots_adjust(bottom=0.5)
    
    # Reverse the colormap to match your get_color inversion
    reversed_cmap = custom_cmap.reversed()
    
    # Now use normal ordering since we reversed the colormap
    norm = plt.Normalize(vmin=-0.2, vmax=0)
    sm = cm.ScalarMappable(cmap=reversed_cmap, norm=norm)
    sm.set_array([])
    
    # Create horizontal colorbar
    cbar = plt.colorbar(sm, cax=ax, orientation='horizontal')
    
    # Labels now directly correspond to your data range
    cbar.set_ticks([-0.2, 0])
    cbar.set_label('Score', labelpad=10)
    
    if show_legend:
        plt.show()
    return fig





def get_color(value): #Gets color for each residue from mean score
    return custom_cmap(1 - value) #Inverts the colors, White is high score and Red is low


def main():
    data = read_scores(file) #Reads data
    region_coords, offset = get_region_info(region)
    filtered, num, coords = pull_scores(data, region_coords) #Gets filtered scores
    residue_values = make_residue_values(filtered, num, coords, offset, analysis) #Makes per-residue mean scores
    normalized_values = normalize_values(residue_values) #Scores normalized to between 0 and 1
    legend = create_colorbar_legend()

    if save_legend:
        legend.savefig('/Users/ivan/Desktop/BARD1_draft_figs/fig5a_BARD1_legend.png', dpi = 500)
    present_residues = list(normalized_values.keys())
    
    for resi in expected_residues:
        if resi not in present_residues:
            print(f'Residue {resi} missing from data, assigning grey color')
            normalized_values[resi] = 1
    #this block does the coloring
    for residue, value in normalized_values.items(): 
        if value == 1:
            color_name = f'color_{chain}_{residue}' #color_A specifies chain A
            cmd.set_color(color_name, [128, 128, 128])  # RGB values
            cmd.color(color_name, f'chain {chain} and resi {residue}') #chain  specifies chain A, change if not chain A
        else:
            color_name = f'color_{chain}_{residue}' #color_A specifies chain A
            color = get_color(value) #Gets color from color map
            cmd.set_color(color_name, [color[0], color[1], color[2]])  # RGB values
            cmd.color(color_name, f'chain {chain} and resi {residue}') #chain  specifies chain A, change if not chain A

    cmd.show('cartoon')

main()