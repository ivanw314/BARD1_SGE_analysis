
"""
This script generates space-filling surfaces in PyMOL for specified structured regions of the BARD1 protein (RING, ARD, BRCT, ARDBRCT) and colors these surfaces based on SGE scores.
Used to generate Fig. 6d and Extended Data Fig. 10 
"""


import pandas as pd
from pymol import cmd
from pymol.cgo import *

sge_scores = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/final_tables/supplementary_file_1_BARD1_SGE_final_table.xlsx' #File of SGE scores
region = 'BRCT' #Structured BARD1 regions to create surface for
chain = 'B' #Specify chain to color 

def region_residues(region): #Takes region input and creates the respective residue numbers
    region_residues = []
    
    if region == 'RING':
        for i in range(26, 123):
            region_residues.append(str(i))
    elif region == 'ARD':
        for i in range(421, 547):
            region_residues.append(str(i))
    elif region == 'BRCT':
        for i in range(566, 778):
            region_residues.append(str(i))
    elif region == 'ARDBRCT':
        for i in range(421, 778):
            region_residues.append(str(i))
    
    return region_residues

def read_scores(file, region_resi): #Reads score file
    df = pd.read_excel(file, sheet_name='scores') #Reads in score file
    
    df = df.loc[df['variant_qc_flag'] != 'WARN'] #Filters out variants with WARN flag
    df = df.rename(columns = {'consequence': 'Consequence', 'amino_acid_change': 'AAsub', 'score': 'snv_score'})
    df = df.loc[df['Consequence'].str.contains('missense_variant')] #Filters only for missense variants

    df['AApos'] = df['AAsub'].transform(lambda x: x[1:-1]) #Creates new amino acid position column

    #Phosphosite site overrides
    phospho_site_override = pd.DataFrame({'pos_id': ['214745805:A','214745805:G', '214745806:A', '214745806:G', '214745806:T' ,'214745805:T', '214745119:A', '214745119:C', '214745119:T', '214745121:A', '214745121:C', '214745121:G', '214745120:T','214745121:A','214745120:A'],
                                           'score': [-0.109977, -0.0175886, -0.0298891, -0.00470975, -0.0413407, 0.0208659, -0.00814722, -0.00228625, 0.00051614, -0.00129351, 0.018338, -0.0955636,  -0.0116846, -0.0049156, -0.0441137],
                                           'amino_acid_change': ['G576V', 'G576A', 'G576C', 'G576R', 'G576S', 'G576D', 'T617T', 'T617T', 'T617T', 'T617S', 'T617A', 'T617P', 'T617N', 'T617S','T617I'],
                                           'functional_consequence': ['functionally_abnormal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'indeterminate', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_abnormal', 'functionally_normal', 'functionally_normal', 'indeterminate'],
                                           'AApos': ['576', '576', '576', '576', '576', '576', '617', '617', '617', '617', '617', '617', '617', '617', '617']
                                          })
    df = pd.concat([df, phospho_site_override]) #Adds phosphosite overrides to score dataframe
    df = df.loc[df['AApos'].isin(region_resi)] #Filters for variants in residues in region of interest
    
    return df

def group_scores(df): #Groups variant scores by AA position and creates calculates min and mean score
    
    df = df[['AApos', 'snv_score']] #Gets AA position and score column only

    grouped = df.groupby('AApos') #Groups by position
    
    min_scores = grouped['snv_score'].min().reset_index() #Summary dataframe with minimum scores
    mean_scores = grouped['snv_score'].mean().reset_index() #Summary dataframe with mean scores
    
    min_scores = min_scores.set_index('AApos')['snv_score'].to_dict() #Turns min scores dataframe into dictionary for coloring
    mean_scores = mean_scores.set_index('AApos')['snv_score'].to_dict() #Turns mean scores dataframe into dictionnary for coloring
  
    return min_scores, mean_scores


def color_surface_by_property(property_dict=None, chain="A", selection="all", palette="rw", 
                            surface_type="surface", transparency=0, surface_name=None): #Creates the colored surfaces
    """
    Colors a molecular surface based on residue properties, creating a separate surface object.
    
    Parameters:
    -----------
    property_dict : dict
        Dictionary mapping residue numbers to property values
    chain : str
        Chain identifier (e.g., "A", "B", etc.)
    selection : str
        Additional PyMOL selection to color (default: "all")
    palette : str
        Color palette to use ("rainbow" or "rw" for red-white)
    surface_type : str
        Type of surface to generate ("surface", "sasurface", "cavitysurface")
    transparency : float
        Surface transparency (0-1, where 1 is fully transparent)
    show_scale : bool
        Whether to show the color scale (default: True)
    surface_name : str
        Name for the new surface object (default: "colored_surface")
    """
    if property_dict is None:
        raise ValueError("Property dictionary must be provided")
    
    if surface_name is None:
        surface_name = "colored_surface"
    
    # Combine chain selection with any additional selection criteria
    full_selection = f"({selection}) and chain {chain}"
    
    # Create a new object for the surface
    cmd.create(surface_name, full_selection)
    
    
    # Show surface only on the new object
    cmd.hide("everything", surface_name)
    cmd.show(surface_type, surface_name)

    # First, color everything gray as default
    cmd.color("gray", surface_name)
    
    # Color mapping function
    def map_value_to_color(value, min_val, max_val, palette):
        if palette == "rainbow":
            spectrum = ["red", "orange", "yellow", "green", "blue", "purple"]
            normalized = (value - min_val) / (max_val - min_val)
            idx = int(normalized * (len(spectrum) - 1))
            return spectrum[idx]
        elif palette in ["rw", "red-white"]:
            # Clamp the value between 0 and 1
            normalized = max(-0.2, min(0, value))
            normalized = 1 + 5 * normalized  #Normalization equation to get normalized variable in [0,1]
            
            # Create a custom color that transitions from red to white
            color_name = f"custom_color_{value}"
            r = 1.0  # Red always stays at 1
            g = b = normalized  # Green and blue increase together from 0 to 1
            
            '''
            # Debug print statements
            if value < 0:
                print(f"Value {value} (< 0): RGB = [{r}, {g}, {b}] - should be pure red [1, 0, 0]")
            elif value > 1:
                print(f"Value {value} (> 1): RGB = [{r}, {g}, {b}] - should be pure white [1, 1, 1]")
            else:
                print(f"Value {value} (between 0-1): RGB = [{r}, {g}, {b}]")
            
            '''
            cmd.set_color(color_name, [r, g, b])

            return color_name
        else:
            raise ValueError("Unsupported palette: choose 'rainbow' or 'rw'")
    
    # Generate color spectrum based on property values
    values = list(property_dict.values())
    min_val, max_val = min(values), max(values)
    
    # Apply colors to the surface object
    for resid, value in property_dict.items():
        color = map_value_to_color(value, min_val, max_val, palette)
        cmd.color(color, f"{surface_name} and resi {resid}")
    
    print(f"Created colored surface object: {surface_name}")
    

def main(chain):
    region_resi = region_residues(region) #Gets region residues
    raw_scores = read_scores(sge_scores, region_resi) #Gets raw SGE scores
    min_scores, mean_scores = group_scores(raw_scores) #Gets min/mean score dataframes

    
    cmd.extend("color_surface_by_property", color_surface_by_property) #Creates color_surface PyMOL command
    
    #color_surface_by_property(property_dict = mean_scores, palette = 'rw', show_scale = True)
    scores = [(min_scores,'SGE_surface_min'), (mean_scores,'SGE_surface_mean')] #Creates an interable list of tuples for min and mean score
    
    for elem in scores: #iterates through the min and mean scores to create both min and mean surfaces
        scores, surface_name = elem
        surface_name = surface_name + "_" + region
        color_surface_by_property(chain = chain, property_dict = scores, palette = 'rw', 
                                       transparency = 0.5, surface_name = surface_name)
main(chain)