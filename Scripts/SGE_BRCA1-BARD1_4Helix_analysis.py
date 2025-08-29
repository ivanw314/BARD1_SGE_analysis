
import pandas as pd
from pymol import cmd

bard1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20250825_BARD1snvscores_filtered.xlsx'
bard1_thresholds = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20250813_BARD1_thresholds.tsv'
brca1_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/20240830_BRCA1_SGE_AllScores.xlsx'


#Get BARD1 thresholds
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

brca1_cutoffs = [-1.328,-0.748]

type = 'min_NP'  # 'min', 'mean', 'min_NP', 'mean_NP' for minimum, mean score, or minimum/mean (proline substituions removed)
chain_info = {'BARD1': 'B', 'BRCA1': 'A'}   #Specify chain identifiers for BARD1 and BRCA1

def read_process_data(bard1_file, brca1_file, type, chains=chain_info):
    bard1_data = pd.read_excel(bard1_file) #Reads BARD1 data from Excel file
    brca1_data = pd.read_excel(brca1_file) #Reads BRCA1 data from Excel file

    #Helix residues for BARD1: [34, 48] and [98, 117]
    #Helix residues for BRCA1: [7, 22] and [80, 97]
    bard1_helix_residues = list(range(26, 48)) + list(range(98, 123)) #BARD1 helix residues
    brca1_helix_residues = list(range(7, 23)) + list(range(80, 98)) #BRCA1 helix residues
 
    bard1_data = bard1_data.rename(columns = {'consequence': 'Consequence'}) #Renames column for consistency
    brca1_data = brca1_data.rename(columns = {'snv_score_minmax': 'score'}) #Renames column for consistency
    bard1_data = bard1_data.loc[bard1_data['Consequence'].isin(['missense_variant'])] #Pulls only missense variants for BARD1
    brca1_data = brca1_data.loc[brca1_data['Consequence'].isin(['missense_variant'])] #pulls only missense variants for BRCA1

    bard1_data['AApos'] = bard1_data['amino_acid_change'].transform(lambda x: x[1:-1]) #Gets amino acid position from the amino acid change

    brca1_data['AApos'] = brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][3:-3]) #Gets amino acid position from the HGVS protein change
    bard1_data = bard1_data.loc[bard1_data['AApos'].astype(int).isin(bard1_helix_residues)] #pulls helix residues for BARD1
    brca1_data = brca1_data.loc[brca1_data['AApos'].astype(int).isin(brca1_helix_residues)] #pulls helix residues for BRCA1

    bard1_data['AAog'] = bard1_data['amino_acid_change'].transform(lambda x: x[0]) #Gets original amino acid for BARD1
    brca1_data['AAog'] = brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][0:3]) #Gets original amino acid for BRCA1
    bard1_data['AAsub'] = bard1_data['amino_acid_change'].transform(lambda x: x[-1]) #Gets amino acid substitution for BARD1
    brca1_data['AAsub'] = brca1_data['hgvs_pro'].transform(lambda x: x[-3:]) #Gets amino acid substitution for BRCA1

    ''' Block to extract data from helix residues and save to Excel files
    brca1_data['functional_consequence'] = 'indeterminate'
    brca1_data.loc[brca1_data['score'] <= brca1_cutoffs[0], 'functional_consequence'] = 'functionally_abnormal'
    brca1_data.loc[brca1_data['score'] >= brca1_cutoffs[1], 'functional_consequence'] = 'functionally_normal'

    genes = [('BARD1', bard1_data), ('BRCA1', brca1_data)]  #List of genes to be used for filtering
    for gene, data in genes:
        file_name = '/Users/ivan/Desktop/test_excel_outputs/20250722_' + gene + '_helix_residues.xlsx'
        data.to_excel(file_name, index=False)

    '''
    #Aggregate data based on the specified type
    if type == 'min':
        bard1_collapsed = bard1_data.groupby('AApos').agg({'score': 'min'}).reset_index()
        brca1_collapsed = brca1_data.groupby('AApos').agg({'score': 'min'}).reset_index()
 
    elif type == 'mean':
        bard1_collapsed = bard1_data.groupby('AApos').agg({'score': 'mean'}).reset_index()
        brca1_collapsed = brca1_data.groupby('AApos').agg({'score': 'mean'}).reset_index()

    elif type == 'min_NP':
        bard1_data = bard1_data.loc[bard1_data['AAsub'] != 'P']
        brca1_data = brca1_data.loc[brca1_data['AAsub'] != 'Pro']

        bard1_collapsed = bard1_data.groupby('AApos').agg({'score': 'min'}).reset_index()
        brca1_collapsed = brca1_data.groupby('AApos').agg({'score': 'min'}).reset_index()

    elif type == 'mean_NP':
        bard1_data = bard1_data.loc[bard1_data['AAsub'] != 'P']
        brca1_data = brca1_data.loc[brca1_data['AAsub'] != 'Pro']

        bard1_collapsed = bard1_data.groupby('AApos').agg({'score': 'mean'}).reset_index()
        brca1_collapsed = brca1_data.groupby('AApos').agg({'score': 'mean'}).reset_index()

    #pulls indeterminate, normal, and abnormal residues based on the cutoffs
    bard1_collapsed['median_consequence'] = 'indeterminate'
    bard1_collapsed.loc[bard1_collapsed['score'] <= bard1_cutoffs[0], 'median_consequence'] = 'functionally_abnormal'
    bard1_collapsed.loc[bard1_collapsed['score'] >= bard1_cutoffs[1], 'median_consequence'] = 'functionally_normal'

    brca1_collapsed['median_consequence'] = 'indeterminate'
    brca1_collapsed.loc[brca1_collapsed['score'] <= brca1_cutoffs[0], 'median_consequence'] = 'functionally_abnormal'
    brca1_collapsed.loc[brca1_collapsed['score'] >= brca1_cutoffs[1], 'median_consequence'] = 'functionally_normal'

    bard1_abnormal = bard1_collapsed.loc[bard1_collapsed['median_consequence'] == 'functionally_abnormal']
    brca1_abnormal = brca1_collapsed.loc[brca1_collapsed['median_consequence'] == 'functionally_abnormal']

    bard1_normal = bard1_collapsed.loc[bard1_collapsed['median_consequence'] == 'functionally_normal']
    brca1_normal = brca1_collapsed.loc[brca1_collapsed['median_consequence'] == 'functionally_normal']

    bard1_indeterminate = bard1_collapsed.loc[bard1_collapsed['median_consequence'] == 'indeterminate']
    brca1_indeterminate = brca1_collapsed.loc[brca1_collapsed['median_consequence'] == 'indeterminate']

    bard1_abnormal_residues = bard1_abnormal['AApos'].astype(int).tolist()
    brca1_abnormal_residues = brca1_abnormal['AApos'].astype(int).tolist()

    bard1_normal_residues = bard1_normal['AApos'].astype(int).tolist()
    brca1_normal_residues = brca1_normal['AApos'].astype(int).tolist()

    bard1_indeterminate_residues = bard1_indeterminate['AApos'].astype(int).tolist()
    brca1_indeterminate_residues = brca1_indeterminate['AApos'].astype(int).tolist()

    to_return = {'BARD1_abnormal': (bard1_abnormal_residues, chain_info['BARD1']),
                 'BARD1_normal': (bard1_normal_residues, chain_info['BARD1']),
                 'BRCA1_abnormal': (brca1_abnormal_residues, chain_info['BRCA1']),
                 'BRCA1_normal': (brca1_normal_residues, chain_info['BRCA1']),
                 'BARD1_indeterminate': (bard1_indeterminate_residues, chain_info['BARD1']),
                 'BRCA1_indeterminate': (brca1_indeterminate_residues, chain_info['BRCA1']),
               } #Final dictionary with residues categorized by their consequence

    return to_return

def generate_selection_string(residues, selection="(all)", chain = chain_info): #Generate selection strings for residues
    """
    Generate a PyMOL selection string for the specified residues.
    Parameters:
    residues (dict): Dictionary with keys containing 'normal'/'abnormal' 
                     and values as (resi_list, chain) tuples.
    selection (str): Base selection to apply residue selection to.
    """
    # Hide all sticks and resets color first
    cmd.hide("sticks", selection)
    all_chains = list(chain_info.values())
    chain_str = "+".join(all_chains)  # Results in "B+A"
    cmd.color('green', f'{selection} and chain {chain_str}')
    
    for key in residues.keys():
        resi_list, chain = residues[key]
        
        if len(resi_list) == 0: #Added to handle cases where no residues are found
            continue
        # Create residue selection string
        resi_str = '+'.join(str(res) for res in resi_list)
        
        # Determine color based on key
        if 'abnormal' in key:
            color = 'red'
        elif 'normal' in key:
            color = 'cyan'
        else:
            color = 'yellow'  # For indeterminate residues
        
        # Create selection name for clarity
        sel_name = f"sel_{key}_{chain}"
        
        # Create the selection
        sel_string = f"({selection}) and chain {chain} and resi {resi_str}"
        
        # Create named selection (helps with debugging)
        cmd.select(sel_name, sel_string)
        
        # Show sticks and color
        cmd.show("sticks", sel_name)
        cmd.color(color, sel_name)
        
        # Optional: print for debugging
        '''
        print(f"Selection '{sel_name}': {sel_string}")
        print(f"Number of atoms selected: {cmd.count_atoms(sel_name)}")
        '''

def main():
    helix_residues = read_process_data(bard1_file, brca1_file, type)

    generate_selection_string(helix_residues)

main()