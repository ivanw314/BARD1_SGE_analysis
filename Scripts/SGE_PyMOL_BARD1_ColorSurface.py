#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:06:12 2025

@author: ivan
"""


import pandas as pd
from pymol import cmd
from pymol.cgo import *

sge_scores = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20250825_BARD1snvscores_filtered.xlsx' #File of SGE scores
region = 'RING' #Structured BARD1 regions to create surface for
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
    df = pd.read_excel(file)
    
    df = df.rename(columns = {'consequence': 'Consequence', 'amino_acid_change': 'AAsub', 'score': 'snv_score'})
    df = df.loc[df['Consequence'].str.contains('missense_variant')] #Filters only for missense variants
    
    df['AApos'] = df['AAsub'].transform(lambda x: x[1:-1]) #Creates new amino acid position column 
    
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
    
    # Color mapping function
    def map_value_to_color(value, min_val, max_val, palette):
        if palette == "rainbow":
            spectrum = ["red", "orange", "yellow", "green", "blue", "purple"]
            normalized = (value - min_val) / (max_val - min_val)
            idx = int(normalized * (len(spectrum) - 1))
            return spectrum[idx]
        elif palette in ["rw", "red-white"]:
            # Clamp the value between 0 and 1
            normalized = max(-0.5, min(0, value))
            normalized = 1 + 2 * normalized  #Normalization equation to get normalized variable in [0,1]
            
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