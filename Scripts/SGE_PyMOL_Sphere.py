import pandas as pd
from pymol import cmd, stored


# === Parameters ===
#Script adapted from S. Fayer
#This script generates figure 6C. Similar figures with spheres can also be generated for all other BARD1 domains and BRCA1 RING/BRCT domains.


base_sphere_size = 0.5  # Starting size for sphere
size_increment = 0.1  # Size increase per variant below threshold
keep_structures = False # Whether to keep existing structures in PyMOL session

# Figure definitions
gene = 'BARD1' #Change to BRCA1 for RING and BRCT comparisons
figure_type = 'ARD' #RING E2 (RING), ARD (ARD H4 interaction), BRCT (BRCT phosphopeptide interaction)

obj_name = f"{gene}_{figure_type}_spheres" #Object name

#Logic for figure types
if figure_type == 'RING':
    print("Generating RING figure")

    if gene == 'BRCA1':
        chain = 'A'
        pdb_id = '1JM7'
        zn_resi = [24, 27, 44, 47, 39, 41, 61, 64]

        sphere_resi = list(range(23, 31)) + list(range(51, 69))
        sphere_resi = [res for res in sphere_resi if res not in zn_resi]  # Exclude zinc-coordinating residues

    elif gene == 'BARD1':
        chain = 'B'
        pdb_id = '1JM7'
        zn_resi = [50, 53, 66, 68, 71, 74, 83, 86]

        sphere_resi = list(range(49, 57)) + list(range(77, 91))
        sphere_resi = [res for res in sphere_resi if res not in zn_resi]  # Exclude zinc-coordinating residues
elif figure_type == 'ARD':
    print("Generating ARD figure")
    pdb_id = '7LYC'
    chain = 'N'
    sphere_resi = list(range(421, 547))

elif figure_type == 'BRCT':
    print("Generating BRCT figure")
    if gene == 'BARD1':
        chain = 'B'
        pdb_id = '3FA2'
        sphere_resi = [575, 576, 617, 619]
    elif gene == 'BRCA1':
        chain = 'A'
        pdb_id = '1T29'
        sphere_resi = [1655, 1656, 1700, 1702]


# Load your score data
if gene == 'BARD1':
    score_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/final_tables/BARD1_SGE_final_table.xlsx'
    score_df = pd.read_excel(score_file, sheet_name='scores')
    score_df = score_df.loc[(~score_df['amino_acid_change'].isin(['---'])) & (score_df['var_type'].isin(['snv']))]
    score_df = score_df.loc[score_df['consequence'].isin(['missense_variant'])]
    score_df = score_df.loc[~score_df['variant_qc_flag'].isin(['WARN'])]

    phospho_site_override = pd.DataFrame({'pos_id': ['214745805:A','214745805:G', '214745806:A', '214745806:G', '214745806:T' ,'214745805:T', '214745119:A', '214745119:C', '214745119:T', '214745121:A', '214745121:C', '214745121:G', '214745120:T','214745121:A','214745120:A'],
                                           'score': [-0.109977, -0.0175886, -0.0298891, -0.00470975, -0.0413407, 0.0208659, -0.00814722, -0.00228625, 0.00051614, -0.00129351, 0.018338, -0.0955636,  -0.0116846, -0.0049156, -0.0441137],
                                           'amino_acid_change': ['G576V', 'G576A', 'G576C', 'G576R', 'G576S', 'G576D', 'T617T', 'T617T', 'T617T', 'T617S', 'T617A', 'T617P', 'T617N', 'T617S','T617I'],
                                           'functional_consequence': ['functionally_abnormal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'indeterminate', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_normal', 'functionally_abnormal', 'functionally_normal', 'functionally_normal', 'indeterminate'],
                                           'AApos': [576, 576, 576, 576, 576, 576, 617, 617, 617, 617, 617, 617, 617, 617, 617]
                                          })

    score_df = pd.concat([score_df, phospho_site_override])
    score_df['aa_pos'] = score_df['amino_acid_change'].transform(lambda x: int(x[1:-1]))
elif gene == 'BRCA1':
    score_file = '/Users/ivan/Documents/GitHub/BARD1_SGE_analysis/Data/extra_data/BRCA1_SGE_data.xlsx'
    old_brca1_data = pd.read_excel(score_file, sheet_name='findlay_2018')
    brca1_data = pd.read_excel(score_file, sheet_name= 'dace_2025')

    brca1_data = brca1_data.rename(columns = {'snv_score_minmax': 'score', 'protPos': 'AApos'}) #Renaming columns for consistency

    old_brca1_min_cutoff = -1.328

    brca1_data = brca1_data.loc[brca1_data['Consequence'].isin(['Missense'])] #pulling only missense variants
    brca1_data['functional_consequence'] = 'indeterminate'  #Setting default function class to indeterminate

    brca1_data.loc[brca1_data['function_class'] == 'LoF', 'functional_consequence'] = 'functionally_abnormal'
    brca1_data.loc[brca1_data['function_class'] == 'Neutral', 'functional_consequence'] = 'functionally_normal'

    old_brca1_data = old_brca1_data.rename(columns = {'snv_score_minmax': 'score'}) #Renaming columns for consistency
    old_brca1_data = old_brca1_data.loc[old_brca1_data['Consequence'].isin(['missense_variant'])] #Pulls only missense variants
    old_brca1_data['AApos'] = old_brca1_data['hgvs_pro'].transform(lambda x: x.split(':')[1].split('.')[1][3:-3]) #Gets amino acid position from the HGVS protein notation
    old_brca1_data['AApos'] = old_brca1_data['AApos'].astype(int) #Ensures amino acid position is an integer
    old_brca1_data['functional_consequence'] = 'indeterminate' #Sets default function class to Intermediate
    old_brca1_data.loc[old_brca1_data['score'] <= old_brca1_min_cutoff, 'functional_consequence'] = 'functionally_abnormal' #Classifies functionally abnormal variants
    old_brca1_data.loc[old_brca1_data['score'] > old_brca1_min_cutoff, 'functional_consequence'] = 'functionally_normal' #Classifies functionally normal variants
    


    new_brca1_pos = list(set(brca1_data['AApos'].tolist())) #List of amino acid positions in the new BRCA1 data
 
    old_brca1_data = old_brca1_data.loc[~(old_brca1_data['AApos'].isin(new_brca1_pos))] #Removes positions already present in the new BRCA1 data


    score_df = pd.concat([brca1_data, old_brca1_data])    #Combines the new and old BRCA1 data
    score_df['AApos'] = score_df['AApos'].astype(int)  #Ensures amino acid position is an integer
    score_df = score_df.rename(columns={'AApos': 'aa_pos'})  #Renames column for consistency
    print(score_df)


score_df = score_df[score_df['aa_pos'].isin(sphere_resi)]
score_df['nAA'] = score_df['amino_acid_change'].transform(lambda x: x[-1])

score_df['LoF_pro'] = 'False'
score_df.loc[(score_df['functional_consequence'] == 'functionally_abnormal') & (score_df['nAA'] == 'P'), 'LoF_pro'] = 'True'

score_df['LoF_pro'] = pd.Categorical(score_df['LoF_pro'], categories=['True', 'False'], ordered=True)

classification_df = score_df.groupby(['aa_pos', 'functional_consequence']).agg(
    count=('LoF_pro', 'size'),
    min_LoF_pro=('LoF_pro', 'min')
).reset_index()

classification_df = classification_df.loc[classification_df['functional_consequence'].isin(['functionally_abnormal'])]

classification_df = classification_df.rename(columns={'aa_pos': 'resi', 'functional_consequence': 'classification', 'count': 'n_variants'})

classification_df['sphere_size'] = classification_df['n_variants'].apply(lambda x: base_sphere_size + (x - 1) * size_increment)

classification_df['color'] = 'gray' # Default color

if figure_type == 'ARD':
    classification_df = classification_df.loc[(classification_df['min_LoF_pro'].isin(['False'])) & (classification_df['n_variants']  >= 2)] #Filtering to only positions with at least one functionally abnormal variant

print(classification_df)

# ============================================================================
# PYMOL VISUALIZATION
# ============================================================================

# Initialize PyMOL
if keep_structures == False:
    cmd.reinitialize()
    cmd.hide("everything")

# Fetch structure
cmd.set('fetch_path', '')  # Don't save fetched files to disk
cmd.fetch(pdb_id)
cmd.create(obj_name, f"{pdb_id} and chain {chain}")
cmd.delete(pdb_id)

# Hide everything first


# Show full chain as cartoon
cmd.show("cartoon", obj_name)
cmd.set("cartoon_fancy_helices", 1)
cmd.color("lightblue", obj_name)

# Add spheres for each position in your dataframe
for idx, row in classification_df.iterrows():
    pos = int(row['resi'])
    sphere_size = float(row['sphere_size'])
    color = row['color']
    # Check if glycine (use CA) or other amino acids (use CB)
    stored.resn_list = []
    cmd.iterate(f"{obj_name} and resi {pos} and name CA", "stored.resn_list.append(resn)")

    if stored.resn_list and stored.resn_list[0] == 'GLY':
        atom_sel = f"{obj_name} and resi {pos} and name CA"
    else:
        atom_sel = f"{obj_name} and resi {pos} and name CB"
    
    # Show sphere with specified size
    cmd.show("spheres", atom_sel)
    cmd.color(color, atom_sel)
    cmd.set("sphere_scale", sphere_size, atom_sel)

# Final styling
cmd.bg_color("white")
cmd.set("ray_shadows", 0)
cmd.set("ambient", 0.4)

print(f"Visualized {len(classification_df)} positions on {pdb_id} chain {chain}")
