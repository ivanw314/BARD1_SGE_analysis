#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:06:12 2025

@author: ivan
"""

import pandas as pd
from pymol import cmd
from pymol.cgo import *


#This script isn't super user-friendly but here's the run down:
#First, please make sure the correct file is selected
#Second, make sure you have the correct region selected, regions are dataset  
#specific and commented out to describe which domain and PDB ID it is from
#Third, 
#BE SURE TO SET THE OFFSET IN PYMOL and set coloring based on chains in the
#make_residue_values() function

file = '/Users/ivan/Desktop/20240830_BRCA1_SGE_AllScores.xlsx' #BRCA1 SGE scores file

region = 'BRCT'
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
        mean = data_filt['snv_score_minmax'].mean() #mean of scores is taken

        #offset needed to assign correct residue number in PyMOL structure
        offset = 1649 #offset in PyMOL structure (what is the number of the AA that starts the coloring?) (1 for RING, 1649 for BRCT)
        
        codon_score[j+ offset] = mean #mean score found and assigned to dictionary with key of base codon number + offset
        
        #this code can be used as a way to use medians for color mapping instead of mean
        #median = data_filt['snv_score_minmax'].median()
        #codon_score[j+1] = median
        j += 1
        
    return codon_score


def color_surface_by_property(property_dict=None, chain="A", selection="all", palette="rwb", 
                            surface_type="surface", transparency=0, show_scale = True):
    """
    Colors a molecular surface based on residue properties.
    
    Parameters:
    -----------
    property_dict : dict
        Dictionary mapping residue numbers to property values
    chain : str
        Chain identifier (e.g., "A", "B", etc.)
    selection : str
        Additional PyMOL selection to color (default: "all")
    palette : str
        Color palette to use ("rainbow", "bwr", "rwb")
    surface_type : str
        Type of surface to generate ("surface", "sasurface", "cavitysurface")
    transparency : float
        Surface transparency (0-1, where 1 is fully transparent)
    """
    if property_dict is None:
        raise ValueError("Property dictionary must be provided")
    
    # Combine chain selection with any additional selection criteria
    full_selection = f"({selection}) and chain {chain}"
    
    # Set up surface settings
    cmd.set("transparency", 0.5)
    cmd.set("surface_quality", 1)  # Ensure good surface quality
    
    # Create the surface for the selected chain
    cmd.show(surface_type, full_selection)
    
    # Set transparency
    cmd.set("transparency", transparency, full_selection)
    
    # Generate color spectrum based on property values
    values = list(property_dict.values())
    min_val, max_val = min(values), max(values)
    
    # Color mapping function
    def map_value_to_color(value, min_val, max_val, palette):
        if palette == "rainbow":
            spectrum = ["red", "orange", "yellow", "green", "blue", "purple"]
            normalized = (value - min_val) / (max_val - min_val)
            idx = int(normalized * (len(spectrum) - 1))
            return spectrum[idx]
        elif palette in ["rw", "red-white"]:
            normalized = (value - min_val) / (max_val - min_val)
            # Create a custom color that transitions from red to white
            color_name = f"custom_color_{value}"
            r = 1.0  # Red always stays at 1
            g = b = normalized  # Green and blue increase together from 0 to 1
            cmd.set_color(color_name, [r, g, b])
            return color_name
        else:
            raise ValueError("Unsupported palette: choose 'rainbow' or 'rw'")
    
    # Apply colors to the selected chain
    for resid, value in property_dict.items():
        color = map_value_to_color(value, min_val, max_val, palette)
        cmd.color(color, f"{full_selection} and resi {resid}")

def color_surface_by_property_beta(property_dict=None, chain="A", selection="all", palette="rw", 
                            surface_type="surface", transparency=0, surface_name=None):
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
            normalized = max(-1.5, min(0, value))
            normalized = 1 + (1/1.5) * normalized  #Normalization equation to get normalized variable in [0,1]
            
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
    

def main():
    data = read_scores(file) #Reads data
    filtered, num, coords = pull_scores(data, regions) #Gets filtered scores
    residue_values = make_residue_values(filtered, num, coords ) #Makes per-residue mean scores
    

    
    cmd.extend("color_surface_by_property", color_surface_by_property)
    
    #color_surface_by_property(property_dict = mean_scores, palette = 'rw', show_scale = True)
    scores = [(residue_values,'SGE_BRCA1_surface_mean')]
    
    for elem in scores:
        scores, surface_name = elem
        surface_name = surface_name + "_" + region
        color_surface_by_property_beta(chain = 'A', property_dict = scores, palette = 'rw', 
                                       transparency = 0.5, surface_name = surface_name)
main()