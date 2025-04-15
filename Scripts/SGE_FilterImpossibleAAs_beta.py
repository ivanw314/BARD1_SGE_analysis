#!/usr/bin/env python
# coding: utf-8

# NEEDS SUPPORT FOR POSITIVE SENSE GENES

# In[9]:


import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq


# In[274]:


oligo_sheet = '/Users/ivan/Desktop/SGE_QC_Test_Files/20240820_CURRENT_SGEoligos_2023.xlsx' #path to SGE oligos 
ref_path = '/Users/ivan/Desktop/SGE_QC_Test_Files/pos_sense_dev/20250415_CTCF_X3_5_ref.xlsx' #path to annotated reference file (you will need to manually give coordinates for where exons are). I originally downloaded the sequence from benchling
sge_scores = '/Users/ivan/Desktop/SGE_QC_Test_Files/pos_sense_dev/CTCF_X3_5.snvscores.tsv'  #path to SGE datafile
gene = 'CTCF' #name of your gene :)
ref_sense = 1 #Sense of the reference file you provide
filtered_file_name = '/Users/ivan/Desktop/20250122_BARD1_PillarProjectScores_filtered' + '.xlsx' #name of saved file


# These first few functions are used to read files or used within other functions 

# In[50]:


def read_unfiltered(sge_file): #read raw SGE data from excel
    df = pd.read_csv(sge_file, sep = '\t')
    df['pos'] = df['pos'].astype(str)
    df['pos_id'] = df['pos'] + ':' + df['allele']

    return df


# In[64]:


def read_ref(ref_file): #reads provided reference file
    ref = pd.read_excel(ref_file)

    return ref


# In[120]:


def read_oligo_sheet(file, gene): #reads the oligo spreadsheet and pulls sheet for specific gene
    oligos_df = pd.read_excel(file, sheet_name = gene)
    original_cols = oligos_df.columns.tolist()
    
    cleaned_cols = []

    i = 0
    while i < 3:
        cleaned_cols.append(original_cols[i].upper())

        i += 1
        
    oligos_df = oligos_df.rename(columns = {original_cols[0]: cleaned_cols[0],
                                            original_cols[1]: cleaned_cols[1],
                                            original_cols[2]: cleaned_cols[2]
                                           }
                                )
    
    oligos_df = oligos_df[['REGION','OLIGO NAME','PRIMER SEQUENCE']]
    
    return oligos_df


# In[20]:


def reverse_complement_string(seq_string): #Reverse complement and returns string
    reverse_seq = seq_string[::-1]
    reverse_comp_list = []
    for char in reverse_seq:
        if char == "A":
            reverse_comp_list.append("T")
        elif char == "G":
            reverse_comp_list.append("C")
        elif char == "C":
            reverse_comp_list.append("G")
        else:
            reverse_comp_list.append("A")
    reverse_compliment_str = "".join(reverse_comp_list)
    return reverse_compliment_str


# In[22]:


def reverse_complement(seq_list): #Reverse complement and returns a list (retains upper/lowercase used in SGE oligo)
    
    i = 0
    #while i < len(seq_list):
        #seq_list[i] = seq_list[i].upper()
        #i += 1
    reverse_seq = seq_list[::-1]
    reverse_comp_list = []
    for char in reverse_seq:
        if char == "A":
            reverse_comp_list.append("T")
        elif char == 'T':
            reverse_comp_list.append("A")
        elif char == "G":
            reverse_comp_list.append("C")
        elif char == "C":
            reverse_comp_list.append("G")
        elif char == 'a':
            reverse_comp_list.append("t")
        elif char == 'c':
            reverse_comp_list.append('g')
        elif char == 'g':
            reverse_comp_list.append('c')
        else:
            reverse_comp_list.append('a')
            
    reverse_compliment_str = "".join(reverse_comp_list)
    return reverse_compliment_str


# In[24]:


def mutate_snvs(dna_sequence): #Mutates all possible SNVs of provided DNA sequence
    snvs = []
    i = 0
    while i < len(dna_sequence):
        if dna_sequence[i] == "A":
            snvs.append(dna_sequence[:i] + "T" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "C" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "G" + dna_sequence[i + 1 :])
        elif dna_sequence[i] == "T":
            snvs.append(dna_sequence[:i] + "A" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "C" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "G" + dna_sequence[i + 1 :])
        elif dna_sequence[i] == "C":
            snvs.append(dna_sequence[:i] + "A" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "T" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "G" + dna_sequence[i + 1 :])
        else:
            snvs.append(dna_sequence[:i] + "A" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "T" + dna_sequence[i + 1 :])
            snvs.append(dna_sequence[:i] + "C" + dna_sequence[i + 1 :])
        i += 1
    return snvs


# These functions below are used to process the oligos spreadsheet, get fixed edit locations, and get data needed for filtering 

# In[215]:


def get_sgeoligos(df, gene, sense): #pulls SGE oligos from spreadsheet and makes reference seq like file
    grouped = df.groupby(by = 'REGION') #grouops by region

    if sense == 0:
        sge_oligos = []
        for target, group in grouped:
            if len(group) > 0:
                split_target = target.split(' ') #rewrites names of targets in oligo spreadsheet to be equal to that in SGE datafile
                region = split_target[1]
                target = gene + '_X'+ region
    
                oligo_name = target + '_sgeoligo'
    
                final_name = gene + '_X' + region.upper() #name equal to name in SGE datafile
                
    
                #starts process to pull coordinates for SGE library
                oligo_coords = group.loc[group['OLIGO NAME'].isin(['SNV library'])]
                coords = str(oligo_coords['PRIMER SEQUENCE'])
                row_split = coords.split(' ')
                row_split = row_split[4]
                coord_split = row_split.split('-')
                chr_coord = coord_split[0]
                
                if chr_coord != 'dtype:': #gets SGE oligo start and end coordinates
                    start_split = chr_coord.split(':')
                    start = int(start_split[1])
                    coord_name = coord_split[1]
                    end_split = coord_name.split('\n')
                    end = int(end_split[0])
                    
                    oligo = group.loc[group['OLIGO NAME'].isin([oligo_name])] #gets SGE oligo
                    oligo = oligo['PRIMER SEQUENCE'].tolist() #puts SGE oligo into list
    
    
                    if len(oligo) > 0:
                        oligo = oligo[0] #gets oligo string
                        oligo = reverse_complement(oligo)
                        coords = [] # list to hold coordinates
                        region = [] # list to hold region name
                        
                        for i in range(start, end + 1): #creates list of coordinates and appends region name each time
                            coords.append(i)
                            region.append(final_name)
                        #coords = coords[::-1] #flips coordinates for antisense gene
        
                        oligo_bp = [] #holds each bp in the library
                        
                        i = 0 
                        while i < len(oligo): #iterates through string and appends letter
                            oligo_bp.append(oligo[i])
                            i += 1
                            
                        region_df = pd.DataFrame({'target': region, 'Reference': oligo_bp, 'pos': coords}) #creates region dataframe
                    sge_oligos.append(region_df) #appends single region oligo dataframe to list
                
        final = pd.concat(sge_oligos) #concatenates all SGE oligo dataframes

    elif sense == 1:
        sge_oligos = []
        for target, group in grouped:
            if len(group) > 0:
                split_target = target.split(' ') #rewrites names of targets in oligo spreadsheet to be equal to that in SGE datafile
                region = split_target[1]
                target = gene + '_X'+ region
    
                oligo_name = target + '_sgeoligo'
    
                final_name = gene + '_X' + region.upper() #name equal to name in SGE datafile
                
                #starts process to pull coordinates for SGE library
                oligo_coords = group.loc[group['OLIGO NAME'].isin(['SNV library'])]
                coords = str(oligo_coords['PRIMER SEQUENCE'])
                row_split = coords.split(' ')
                row_split = row_split[4]
                coord_split = row_split.split('-')
                chr_coord = coord_split[0]
                
                if chr_coord != 'dtype:': #gets SGE oligo start and end coordinates
                    start_split = chr_coord.split(':')
                    start = int(start_split[1])
                    coord_name = coord_split[1]
                    end_split = coord_name.split('\n')
                    end = int(end_split[0])
                    
                    oligo = group.loc[group['OLIGO NAME'].isin([oligo_name])] #gets SGE oligo
                    oligo = oligo['PRIMER SEQUENCE'].tolist() #puts SGE oligo into list
    
    
                    if len(oligo) > 0:
                        oligo = oligo[0] #gets oligo string
                        coords = [] # list to hold coordinates
                        region = [] # list to hold region name
                        
                        for i in range(start, end + 1): #creates list of coordinates and appends region name each time
                            coords.append(i)
                            region.append(final_name)
        
                        oligo_bp = [] #holds each bp in the library
                        
                        i = 0 
                        while i < len(oligo): #iterates through string and appends letter
                            oligo_bp.append(oligo[i])
                            i += 1
                        #print(len(region), len(oligo_bp), len(coords)) #QC to check length of arrays
                        region_df = pd.DataFrame({'target': region, 'Reference': oligo_bp, 'pos': coords}) #creates region dataframe
                        
                    sge_oligos.append(region_df) #appends single region oligo dataframe to list
                
        final = pd.concat(sge_oligos) #concatenates all SGE oligo dataframes


    return final


# In[205]:


def get_edits(oligos,gene): #gets fixed edits for that target
    grouped = oligos.groupby(by = 'REGION') #groups dataframe by SGE region
    edit_dicts = []
    
    for target, grouped in grouped:

        split_target = target.split(' ') #rewrites names of targets in oligo spreadsheet to be equal to that in SGE datafile
        region = split_target[1]
        region = region.upper()
        target = gene + '_X'+ region

        edits = grouped.loc[grouped['OLIGO NAME'].isin(['Edits'])] #does string things to extract edits from spreadsheet 
        text_edit = str(edits['PRIMER SEQUENCE'])
        colon_split_edit = text_edit.split(':') #splits on the shared colon
        coord_edits = colon_split_edit[1] #Split string gives list, gets list component with coords of edits
        comma_coord_split = coord_edits.split(',') #Edits are split by comma, splits the edits by the comma
        
        edit_list = [] #list to hold edits for this target
        edit_dict = {} #placeholder dictionary to store target name and edit_list
        for elem in comma_coord_split: #iterates through the list created by splitting the string containing the edits
            chars = [] #list to hold each number
            for char in elem: #loop checks to see if each character is a number
                if char.isdigit():
                    chars.append(char)
                    
            edit_coord = ''.join(chars) #joins list of characters to yield complete corrdinate
            if len(edit_coord) > 1:
                edit_list.append(int(edit_coord)) #appends coordinate to list
                edit_dict[target] = edit_list #creates dictionary of edits for that target
        if len(edit_dict) == 1:
            edit_dicts.append(edit_dict) #appends target-specific dictionary to list containing edits for all regions
    
    return edit_dicts


# In[249]:


def prep_filter_type(dicts,ref_path,sense): #generates dataframe that will be used to determine what kind of fixed edit is present

    ref_all = pd.read_excel(ref_path) #reads reference file
    
    #lists for final dataframe
    target_list = [] #holds target names
    pos_list = [] #holds edit coordinates
    ref_anti = [] #holds the antisense 4bp window
    ref_sense = [] #holds the sense 4bp window
    
    for elem in dicts: 
        keys = elem.keys()
        for key in keys:
            target = key #stores name of SGE region 
            
        ref_target = ref_all.loc[ref_all['target'].isin([target])] #pulls reference sequence for that SGE library  
        ref_target = ref_target.loc[ref_target['Intron/Exon'].isin(['Exon'])] #excludes bases annotated as intronic to remove intronic HDR markers
        
        edits = elem[target] #gets list of edits

        for elem in edits:
                if sense == 0:
                    target_ref_a = ref_target
                    target_ref_s = ref_target
                    min_2 = elem - 2
    
                    edit_coords_0 = []
                    for i in range(min_2 + 1, elem + 3): #used to get 4bp window on antisense strand
                        edit_coords_0.append(i)
                    edit_coords_1 = []
                    for i in range(min_2, elem + 2): #used to get 4 bp window on sense strand
                        edit_coords_1.append(i)
                            
                    ref_0 = target_ref_a.loc[target_ref_a['pos'].isin(edit_coords_0)]
                    ref_0 = ref_0['Reference'].to_list()
                    ref_0 = reverse_complement(ref_0) #Gets 4 bp window that contains 2 bp upstream of edit and 1 bp downstream on antisense strand
    
                    ref_1 = target_ref_s.loc[target_ref_s['pos'].isin(edit_coords_1)]
                    ref_1 = ref_1['Reference'].to_list()
                    ref_1 = ''.join(ref_1) #Gets 4 bp window that contains 2 bp upstream of edit and 1 bp downstream on sense strand


                    if len(ref_0) > 0: #tests to see if any bp are pulled (0 for intronic edits)
                        target_list.append(target)
                        pos_list.append(elem)
                        ref_anti.append(ref_0)
                        ref_sense.append(ref_1)
                    
                elif sense == 1:
                    target_ref_a = ref_target
                    target_ref_s = ref_target
                    min_2 = elem - 2
    
                    edit_coords_0 = []
                    for i in range(min_2 + 1, elem + 3): #used to get 4bp window on antisense strand
                        edit_coords_0.append(i)
                    edit_coords_1 = []
                    for i in range(min_2, elem + 2): #used to get 4 bp window on sense strand
                        edit_coords_1.append(i)
                            
                    ref_0 = target_ref_a.loc[target_ref_a['pos'].isin(edit_coords_0)]
                    ref_0 = ref_0['Reference'].to_list()
                    ref_0 = reverse_complement(ref_0) #Gets 4 bp window that contains 2 bp upstream of edit and 1 bp downstream on antisense strand
    
                    ref_1 = target_ref_s.loc[target_ref_s['pos'].isin(edit_coords_1)]
                    ref_1 = ref_1['Reference'].to_list()
                    ref_1 = ''.join(ref_1) #Gets 4 bp window that contains 2 bp upstream of edit and 1 bp downstream on sense strand


                    if len(ref_0) > 0: #tests to see if any bp are pulled (0 for intronic edits)
                        target_list.append(target)
                        pos_list.append(elem)
                        ref_anti.append(ref_0)
                        ref_sense.append(ref_1)
                    
                

    df = pd.DataFrame({'target': target_list, 'edit_coord': pos_list, 'anti': ref_anti, 'sense': ref_sense})

    return df, ref_all

            


# In[278]:


def type_edits(df,ref, oligos, sense): #determines what type of edit is present in preparation for filtering. Needs prep_filter dataframe, ref seq dataframe, SGE oligos dataframe, and gene sense

    #Lists to store characterized fixed edit tuples in form of (SGE region, coordinate of edit, bp of PAM site impacted, sense of PAM site, base change on sense strand)
    pams = [] #list to store PAM edits (edits impacting NGG sites)
    edits = [] #list to store other fixed edits (not at NGG sites, adjacent fixed edits under guide sequence)
    same_codon_doubles = [] #list to store coordinates containing two fixed edits at the same codon
    
    grouped = df.groupby(by = 'target') #groups the dataframe created by the prep_filter_type function by SGE region
    non_pam_coords  = [] #list to hold coordinates of fixed edits not at NGG sites
    
    for target, group in grouped: #iterates through each region and determines what fixed edits are present 
        target_edits = group['edit_coord'].tolist() #list of fixed edits for that SGE region
        pam_df = group #group is a dataframe and is also assigned to this variable for future use in determining what PAMs are present
        ref_target = ref.loc[ref['target'].isin([target])] #pulls out reference sequence for SGE region
        target_oligo = oligos.loc[oligos['target'].isin([target])] #pulls out specific SGE oligo from concatenated dataframe 
        
        i = 0
        while (i + 1) < len(target_edits): #iterates through list of fixed edits and sorts all edits except for those at NGG PAM sites
            test = target_edits[i + 1] - target_edits[i] #initial test to sort same codon edits and adjancent edits under guide

            if test == 2: #Fixed edits within the same codon have difference between genomic coordinates of 2
                tuple = (target, target_edits[i], target_edits[i + 1]) #creates tuple with (SGE region, 1st fixed edit, 2nd fixed edit)

                #appends coordinates to non PAM site edits list
                non_pam_coords.append(target_edits[i])
                non_pam_coords.append(target_edits[i + 1])

                #appends tuple to holding list for same codon fixed edits
                same_codon_doubles.append(tuple)
            
            elif test == 3: #Adjacent fixed edits under the guide sequence will have difference between genomic coordinates of 3
                edit_one = int(target_edits[i]) #Gets coordinate of first fixed edit
                edit_two = int(target_edits[i + 1]) #Gets coordinate of second fixed edit

                #lists to hold the coordinates of the adjacent codons impacted by fixed edits
                codon_one_coords = []
                codon_two_coords = []

                #for loops used to create the genomic coordinates for each codon
                for j in range(edit_one, (edit_one + 3)):
                    codon_one_coords.append(j)
                for j in range(edit_two, edit_two + 3):
                    codon_two_coords.append(j)

                #Gets the codons impacted by fixed edits from reference and the SGE oligo 
                codon1_df = ref_target.loc[ref_target['pos'].isin(codon_one_coords)]
                codon2_df = ref_target.loc[ref_target['pos'].isin(codon_two_coords)]
                codon1_oligo = target_oligo.loc[target_oligo['pos'].isin(codon_one_coords)]
                codon2_oligo = target_oligo.loc[target_oligo['pos'].isin(codon_two_coords)]

                #For both codons, dataframes are merged
                codon1_merged = pd.merge(codon1_df,codon1_oligo, how = 'outer', on = ['target', 'Reference', 'pos'], indicator = True)
                codon2_merged = pd.merge(codon2_df,codon2_oligo, how = 'outer', on = ['target', 'Reference', 'pos'], indicator = True)

                #Dataframes are filtered to include only the base that is different (i.e. what the fixed edit is)
                codon1_notshared = codon1_merged.loc[~codon1_merged['_merge'].isin(['both'])]
                codon2_notshared = codon2_merged.loc[~codon2_merged['_merge'].isin(['both'])]

                #Concatenates the dataframes for this pair of fixed edits
                concat_codons = pd.concat([codon1_notshared,codon2_notshared])
                concat_codons = concat_codons.drop(columns = ['Unnamed: 0','Intron/Exon', '_merge'])

                #Groups fixed edits by genomic coordinate in preparation for tuple formation
                grouped_codons = concat_codons.groupby( by = 'pos')

                #iterates through group object and creates tuple for filtering 
                for pos, group in grouped_codons:
                    name = group['target'][3]
                    edit = group['Reference'][3]
                    edit = edit.upper()
                    pos = group['pos'][3]
                    tuple = (name, pos,3,0, edit) #PAM position and sense component of this tuple not used and thus are fixed at 3 and 0
                    non_pam_coords.append(pos)
                    edits.append(tuple)


            i += 1

        #From pre-filtering dataframe, retrieves fixed eidts at PAM sites
        pam_df = pam_df.loc[~pam_df['edit_coord'].isin(non_pam_coords)]
        pam_df = pam_df.reset_index(drop = True) #resets index for while loop 

        #Sorts fixed PAM edits
        if len(pam_df) > 0: #if statement needed to filter out empty dataframes generates from fixed edits not at PAM sites
            i = 0
            pattern = r"[ACTG]GG" #Regular expression for recognizing PAM sites

            #iterates through dataframe row by row and sorts PAM site
            while i < len(pam_df):
                pos = pam_df['edit_coord'][i] #genomic coordinate of fixed edit
                anti = (pam_df['anti'][i]).upper() #antisense 4bp window created from pre-filtering function
                sense = (pam_df['sense'][i]).upper() #sense 4bp window created from pre-filtering function

                if anti == 'GGGG' or sense == 'GGGG': #if any 4bp window contains GGGG, then it cannot be automatically filtered, user entry required
                    print('Exception: GGGG PAM in ', target)
                    pos_changed = int(input('PAM position changed: ')) #position in PAM that the fixed edit occurred (2 or 3)
                    pam_sense = int(input('Sense of PAM: ')) #is the PAM on sense (1) or antisense (0) strand

                    #if and elif statements determine what basechange occurs on the sense strand (used in SGE data output) occurs based on sense of PAM site
                    if pam_sense == 0:
                        base_change = 'G'
                    elif pam_sense == 1:
                        base_change = 'C'

                    #creates and appends tuple with PAM information
                    tuple = (target, pos, pos_changed, pam_sense, base_change)
                    pams.append(tuple)

                #Using regex, searches antisense windows for PAMs modified in the 3rd position on antisense strand
                elif re.search(pattern, anti[0:3]):
                    pos_changed = 3
                    pam_sense = 0
                    base_change = 'G' #base change is G as G->C on antisense strand leads to C->G change on sense strand
                    tuple = (target, pos, pos_changed, pam_sense, base_change)
                    pams.append(tuple)

                #Using regex, searches antisense windows for PAMs modified in the 2nd position on antisense strand
                elif re.search(pattern, anti[1:4]):
                    pos_changed = 2
                    pam_sense = 0
                    base_change = 'G'
                    tuple = (target, pos, pos_changed, pam_sense, base_change)
                    pams.append(tuple)

                #Using regex, searches sense windows for PAMs modified in the 3rd position on sense strand
                elif re.search(pattern, sense[0:3]):
                    pos_changed = 3
                    pam_sense = 1
                    base_change = 'C'
                    tuple = (target, pos, pos_changed, pam_sense, base_change)
                    pams.append(tuple)

                #Using regex, searches sense windows for PAMs modified in the 2nd position on sense strand
                elif re.search(pattern, sense[1:4]):
                    pos_changed = 2
                    pam_sense = 1
                    base_change = 'C'
                    tuple = (target, pos, pos_changed, pam_sense, base_change)
                    pams.append(tuple)

                i += 1
                    
 
    return pams, edits, same_codon_doubles
    


# These functions are used for filtering

# In[170]:


def remove_low_freq(df): #Removes all variants with low frequency either in SNV library or at Day 5 (Deprecated 4/14/25)
    
    #remove SNVs that were in library less than 1 in 10,0000
    df_snv = df.loc[df['snvlib_freq'] > 0.001]

    #remove variants present at day 5 with frequency less than 1 in 100,000
    df_freq_1 = df_snv.loc[(df_snv['D05_R1R2R3_freq'] > 0.0001) & (df_snv['D05_R4R5R6_freq'] > 0.0001)]
    #df_freq_2 = df_snv.loc[(df_snv['D05_R1R4_freq'] > 0.0001) & (df_snv['D05_R2R5_freq'] > 0.0001) & (df_snv['D05_R3R6_freq'] > 0.0001)]
    df_freq_3 = df_snv.loc[(df_snv['D05_R1R4_freq'] > 0.0001) & (df_snv['D05_R2R5_freq'] > 0.0001)]

    to_concat = [df_freq_1, df_freq_3]
    df_freq = pd.concat(to_concat)
    
    return df_freq


# In[304]:


def pam_filter(df, pams, sense): #Iterates through fixed edits at NGG sites and removes new PAMs one bp away from original PAM as well as variants at PAM sites

    for elem in pams: #iterates through list of PAM edits
        var_to_filter = []
        target, coord, pos, sense,fixed = elem

        #filters out all variants at PAM edits
        var_to_filter.append(str(coord) + ':G')
        var_to_filter.append(str(coord) + ':C')
        var_to_filter.append(str(coord) + ':A')
        var_to_filter.append(str(coord) + ':T')

        #Filters out new NGG PAMs one bp away from original PAM (possible with position 2 edit)
        if pos == 2:
            if sense == 1: #for PAM on positive sense strand
                new_pam = coord +  2 #plus 2 for antisense gene
                to_filter = str(new_pam) + ':G'
                var_to_filter.append(to_filter)
            if sense == 0: #for PAM on negative sense strand
                new_pam = coord -  2 # -2 for antisense gene
                to_filter = str(new_pam) + ':C' # ':C' as a N>C change on - sense will yield G on + sense strand
                var_to_filter.append(to_filter)
           
        df = df.loc[~((df['pos_id'].isin(var_to_filter)) & (df['target'].isin([target])))] #filters out new PAM and PAM site edits present in the same SGE target

    return df


# In[497]:


def get_impossible_snvs(pam_filtered,all_fixed_edits,ref, gene_sense):
    
    #ref = ref.drop_duplicates(subset = ['pos'], keep = 'last') #drops duplicates from reference
    impossible_snvs = []
    impossible_snvs_df = []
    
    if gene_sense == 0:
        for elem in all_fixed_edits: #gets list of possible amino acid changes
            ref_mod = ref
            target, coord, pos, sense, fixed = elem #unpacks edit tuple provided by user
            ref_mod = ref_mod.loc[ref_mod['target'].isin([target])]
            codon_coords = [] #stores coordinates for codon
            
            for i in range(coord, coord + 3): #creates coordinates for codon
                codon_coords.append(i)

            codon_coords = codon_coords[::-1] #Reverses order of coordinate list for antisense gene
            
            codon_ref = ref_mod.loc[ref_mod['pos'].isin(codon_coords)] #pulls out reference sequence from reference sequence dataframe
            codon = codon_ref['Reference'].tolist() #turns reference sequence to list
            codon = ''.join(codon) #joins references sequence to string
            antisense_codon = reverse_complement_string(codon) #creates coding sequence

            possible_snvs = mutate_snvs(antisense_codon) #generates all possible SNVs

            if len(codon_ref) % 3 != 0 or len(codon_ref) != 3: #error checking - prints out outputs from previous that are not 1 codon only
                print(codon_ref)
                
            #creates and stores all possible AA changes from reference sequence
            possible_vars = []
            for elem in possible_snvs:
                var = Seq(elem)
                aa = var.translate()
                possible_vars.append(str(aa))
                  
            fixed_coords = []
            fixed_vars = []
            for i in range(coord + 1, coord + 3): #gets list of amino acid changes in context of fixed edits
                fixed_coords.append(i)
            rest_codon_fix = ref_mod.loc[ref_mod['pos'].isin(fixed_coords)] #gets bases in codon excluding fixed edit from reference dataframe

            rest_codon_fix = rest_codon_fix['Reference'].tolist() #reference dataframe to list
            codon_fix = [fixed] #fixed edit only
            codon_fix += rest_codon_fix #full codon created
            
            codon_fix = ''.join(codon_fix)
            fixed_coding = reverse_complement_string(codon_fix) #reverse complements for genes on negative sense strand
            fixed_edit = fixed_coding[2] #gets fixed edit 
            to_mutate = fixed_coding[0:2] #gets not fixed edit
            fixed_snvs = mutate_snvs(to_mutate) #creates all SNVs for not fixed edits
            for elem in fixed_snvs: #generates all possible SNVs with the fixed edit
                var = elem + fixed_edit
                fixed_vars.append(var)

            #Creates list of impossible variants (one list per fixed edit)
            impossible = [] 
            j = 0
            while j < len(fixed_vars):
                fixed = Seq(fixed_vars[j]) #gets a variant with fixed edit
                if len(fixed) != 3:
                    print(fixed)
                fixed_aa = str(fixed.translate()) #translates to amino acid
                if fixed_aa in possible_vars:
                    j +=1
                else: #if amino acid isn't in the possible list, it is appended to the impossible list
                    impossible.append(str(fixed))
                    j +=1
                    
                
            #print(impossible)
            
            #takes impossible SNVs and creates position IDs to filter out 
            
            #takes unmutated, fixed edit codon, splits up by basepair
            split_codon = [] #for ref column in final df that goes to list
            k = 0
            while k < len(fixed_coding):
                split_codon.append(str(fixed_coding)[k])
                k += 1

            #iterates through each list of impossible SNVs to create position IDs for each list to store in df
            for elem in impossible:
                coords = codon_coords #coordinates of codon
                ref_codon = split_codon #reference is the WT codon with fixed edit
                fixed_codon = [] #list to hold basepairs of impossible codon
                pos_id_list = [] #list to hold position IDs of impossible codons for this fixed edit

                #splits up impossible codon by bp
                i = 0
                while i < len(elem):
                    fixed_codon.append(elem[i])
                    i += 1

                #iterates through each basepair to determine where variant occurs and gets coordinate to create position ID
                j = 0
                while j < len(fixed_codon):
                    if fixed_codon[j] != ref_codon[j]: #if statement determines mismatch
                        coord = coords[j] #gets coordinate
                        pos_id = str(coord) + ':' + reverse_complement(fixed_codon[j]) #creates position ID (reverse complement needed due to pos_ids being reported on sense strand)
                        pos_id_list.append(pos_id)
                    j += 1

                #creates dictionary then dataframe with position IDs for this fixed edit
                data = {'pos_id': pos_id_list} 
                target_names = []
                for i in range(len(data)):
                    target_names.append(target)
                data['target'] = target_names   
                df = pd.DataFrame(data)
                
                #print(df)
                impossible_snvs.append(df)

        impossible_snvs_df = pd.concat(impossible_snvs) #concatenates all impossible SNV dataframes

    elif gene_sense == 1:
        
        for elem in all_fixed_edits: #gets list of possible amino acid changes
            ref_mod = ref
            target, coord, pos, sense, fixed = elem #unpacks edit tuple provided by user
            ref_mod = ref_mod.loc[ref_mod['target'].isin([target])]
            codon_coords = [] #stores coordinates for codon
            for i in range(coord - 2, coord + 1): #creates coordinates for codon
                codon_coords.append(i)
            
            codon_ref = ref_mod.loc[ref_mod['pos'].isin(codon_coords)] #pulls out reference sequence from reference sequence dataframe
            codon = codon_ref['Reference'].tolist() #turns reference sequence to list
            codon = ''.join(codon) #joins references sequence to string


            possible_snvs = mutate_snvs(codon) #generates all possible SNVs

            if len(codon_ref) % 3 != 0 or len(codon_ref) != 3: #error checking - prints out outputs from previous that are not 1 codon only
                print(codon_ref)
                
            #creates and stores all possible AA changes from reference sequence
            possible_vars = []
            for elem in possible_snvs:
                var = Seq(elem)
                aa = var.translate()
                possible_vars.append(str(aa))
                
            
            fixed_coords = []
            fixed_vars = []
            for i in range(coord - 2, coord): #gets list of amino acid changes in context of fixed edits
                fixed_coords.append(i)
            rest_codon_fix = ref_mod.loc[ref_mod['pos'].isin(fixed_coords)] #gets bases in codon excluding fixed edit from reference dataframe

            rest_codon_fix = rest_codon_fix['Reference'].tolist() #reference dataframe to list
            codon_fix = [fixed] #fixed edit only
            codon_fix = rest_codon_fix + codon_fix #full codon created
            
            fixed_coding = ''.join(codon_fix) #joins codon into one string

            fixed_edit = fixed_coding[2] #gets fixed edit 
            to_mutate = fixed_coding[0:2] #gets not fixed edit

            fixed_snvs = mutate_snvs(to_mutate) #creates all SNVs for not fixed edits
            for elem in fixed_snvs: #generates all possible SNVs with the fixed edit
                var = elem + fixed_edit
                fixed_vars.append(var)
                

            #Creates list of impossible variants (one list per fixed edit)
            impossible = [] 
            j = 0
            while j < len(fixed_vars):
                fixed = Seq(fixed_vars[j]) #gets a variant with fixed edit
                if len(fixed) != 3:
                    print(fixed)
                fixed_aa = str(fixed.translate()) #translates to amino acid
                if fixed_aa in possible_vars:
                    j +=1
                else: #if amino acid isn't in the possible list, it is appended to the impossible list
                    impossible.append(str(fixed))
                    j +=1
                    
                
            #print(impossible)
            
            #takes impossible SNVs and creates position IDs to filter out 
            
            #takes unmutated, fixed edit codon, splits up by basepair
            split_codon = [] #for ref column in final df that goes to list
            k = 0
            while k < len(fixed_coding):
                split_codon.append(str(fixed_coding)[k])
                k += 1

            #iterates through each list of impossible SNVs to create position IDs for each list to store in df
            for elem in impossible:
                coords = codon_coords #coordinates of codon
                ref_codon = split_codon #reference is the WT codon with fixed edit
                fixed_codon = [] #list to hold basepairs of impossible codon
                pos_id_list = [] #list to hold position IDs of impossible codons for this fixed edit

                #splits up impossible codon by bp
                i = 0
                while i < len(elem):
                    fixed_codon.append(elem[i])
                    i += 1

                #iterates through each basepair to determine where variant occurs and gets coordinate to create position ID
                j = 0
                while j < len(fixed_codon):
                    if fixed_codon[j] != ref_codon[j]: #if statement determines mismatch
                        coord = coords[j] #gets coordinate
                        pos_id = str(coord) + ':' + fixed_codon[j] #creates position ID (reverse complement needed due to pos_ids being reported on sense strand)
                        pos_id_list.append(pos_id)
                    j += 1

                #creates dictionary then dataframe with position IDs for this fixed edit
                data = {'pos_id': pos_id_list} 
                target_names = []
                for i in range(len(data)):
                    target_names.append(target)
                data['target'] = target_names   
                df = pd.DataFrame(data)
            
                impossible_snvs.append(df)

        impossible_snvs_df = pd.concat(impossible_snvs) #concatenates all impossible SNV dataframes
        
    return impossible_snvs_df


# In[470]:


def filter_impossible_snvs(data_df, impossible_df): #Removes the impossible missense variants created due to context of fixed edits

    grouped = impossible_df.groupby(by = 'target') #Groups dataframe created by previous function by SGE region

    for target, group in grouped: #iterates through each group
        to_remove = group['pos_id'].tolist() #each position ID to remove in that region is turned into list

        data_df = data_df.loc[~((data_df['pos_id'].isin(to_remove)) & (data_df['target'].isin([target])))] #Datapoints with the same target and position ID are filtered out

    data_filtered = data_df


    return data_filtered


# In[472]:


def remove_doubles(df,doubles): #iterates through list of same codon double edits and filters out all variants within those codons

    for elem in doubles: #iterates through list of same codon doubles
        target, start, end = elem #unpacks tuple
        double_coords = [] #list to hold codons with double edits
        
        for i in range(start, end + 1):
            double_coords.append(i)
        df = df.loc[~((df['pos'].isin(double_coords)) & (df['target'].isin([target])))] #removes double-edit codons
    
    #reindexes
    df = df.reset_index(drop = True)

    return df


# In[492]:


def main(): #Runs all functions
    
    #Pre-filtering functions
    oligos = read_oligo_sheet(oligo_sheet,gene)
    sge_oligos = get_sgeoligos(oligos,gene,ref_sense)
    edits = get_edits(oligos, gene)
    type_ready, reference = prep_filter_type(edits,ref_path,ref_sense)
    pams,edits, doubles = type_edits(type_ready, reference, sge_oligos, ref_sense)
    all_fixed_edits = pams + edits
    #displays fixed edits
    print('Detected Edits follow: ', '\n', 
          'PAM Edits: ',pams, '\n', 
          'Other Fixed Edits: ', edits, '\n',
          'Codon Double Edits: ',doubles)
    
    #Functions for filtering
    raw_data = read_unfiltered(sge_scores)
    ref = read_ref(ref_path)
    pam_filtered = pam_filter(raw_data,pams, ref_sense)
    impossible = get_impossible_snvs(pam_filtered, all_fixed_edits,ref,ref_sense)
    snvs_filtered = filter_impossible_snvs(pam_filtered, impossible)
    filtered_data = remove_doubles(snvs_filtered,doubles)

    #SAVES FILE (Comment out if you don't want it!)
    #filtered_data.to_excel(filtered_file_name, index = False)

    #Displays final filtered data and number of datapoints removed at each step
    #print(filtered_data)
    print('Filtering statistics: ', '\n',
          'RAW: ', len(raw_data), '\n',
          'POST PAM FILTER: ', len(pam_filtered), '\n',
          'POST SNV FILTER: ', len(snvs_filtered),'\n',
          'DONE: ', len(filtered_data))


# In[494]:


main()


# In[ ]:





# In[ ]:




