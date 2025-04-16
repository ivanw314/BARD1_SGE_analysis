# BARD1_SGE_analysis

This repository contains all Python notebooks and scripts used to filter data and create figures used in the BARD1 SGE paper. All notebooks/scripts are organized into one of two folders. Notebooks and Scripts. Datasets used to generate figures are located in the Data folder.

## Notebooks
The Notebooks folder contains Python Notebooks that create individual panels used in the final figures. The vizualization that will be created by each notebook is noted in the name of the notebook. For the intricacies of each notebook, notebooks have been commented out.

These notebooks are:
* SGE_Histogram - Creates histogram of all BARD1 SGE scores **(Needs SGE Data File)**
* SGE_RidegelinePlot_ByVarType_withTicks - Creates ridgeline plot of all scores separated by mutation type with ticks for each individual datapoint **(Needs SGE Data File)**
* SGE_StackedBars_byExon - Stacked bar chart showing fraction functional/non-functional variants across all BARD1 exons facted by missense/stop/synonymous variants **(Needs SGE Data File)**
* SGE_AAHeatMap - Heat map of all amino acid substitutions assayed by SGE **(Needs SGE Data File)**
* SGE_CorrelationMatrix_Heatmap - Correlation heatmap between all replicates **(Needs SGE Data File)**
* SGE_DeletionFigure_ReadDepth - Line plot of median read depth across all coding regions in BARD1 **(Needs: BARD1 Depth file and Deletion Input File)**
* SGE_BARD1_Structure_Faceted - Generates 2 plots. Violin plot showing distribution of scores across the 3 strucutred domains of BARD1. Interactive histograms showing distirbution of scores within each strucutred domain.  **(Needs SGE Data File)**
* SGE_Basic_ClinVar_analysis - Generates 3 plots. 2 histograms showing distribution of current BARD1 variants in ClinVar. 1 Ridgeline plot showing distribution of SGE scores for each ClinVar classification and variants in gnomAD. **(Needs: SGE Data File, BARD1 ClinVar File, and BARD1 gnomAD file)**
* SGE_OrthogonalAssays - Generates 2 plots. Histogram and dotplot showing SGE scores of BARD1 variants assayed in orthogonal assays. **(Needs: SGE Data file and BARD1 Orthogonal Assays File)**
* SGE_gnomAD_MAF_vs_Score - Creates scatterplot of gnomAD mean allele frequency vs. SGE score **(Needs: SGE Data File and BARD1 gnomAD file)**
* SGE_ATGlib_analysis - Generates continuous log-fold change scores for variants in the BARD1_X1A double start codon mutant experiment as well as basic visualizations of the results **(Needs ATG counts folder and X1A annotations file)**
* SGE_ATGlib_HeatMap - Generates heat map visualizations of BARD1_X1A double start codon mutant experiment faceted by canonical start **(Needs X1A Start Score File)**
  
## Scripts
The Scripts folder contains two types of files: Python scripts used to generated figures in PyMOL and the Python notebook used to filter all data prior to figure generation. All scripts are labeled with what visualization they will create once run in PyMOL. Specifics for running each script are included within commented annotations in each individual script. 

These scripts are:
* SGE_BARD1_colorPyMOL_MIS_only - Colors ribbon cartoon protein structures using missense SGE scores  **(Needs SGE Data File)**
* SGE_PyMOL_BARD1_ColorSurface - Generates colored surfaces displaying minimum and mean missense surface score  **(Needs SGE Data File)**
* SGE_PyMOL_BRCA1_ColorSurface - Generates colored surfaces displaying minimum and mean missense surface score. Specific script to handle previous BRCA1 SGE data **(Needs BRCA1 SGE File)**s
* SGE_DataFilter_ReferenceGenerator - Generates library reference file needed for filtering impossible amino acid changes **(Needs BARD1 fasta file and BARD1 input file)**
* SGE_FilterImpossibleAAs_beta - Filters out variants that are likely to yield poor data. Includes: amino acid changes not possible in WT sequence and variants that craete new NGG PAM site near guide site. Also includes functions for handling removing variants that are low-frequency, but not currently in use. **(Needs: SGE Data File, oligo design file, and lib reference file)**


## Data
This folder contains information needed to run notebooks for figure generation and filtering of variants to yield final data set. This includes the raw SGE scores in TSV format as well as data from databases such as ClinVar and gnomAD. All data files dated based on when they were retrieved. 

These datafiles are:
* 20250122_BARD1_SGEscores_wAAsub.xlsx - BARD1 SGE scores based on LFC with column for amino acid substitution. **(SGE Data File)**
* 20240820_CURRENT_SGEoligos_2023.xlsx - All designed SGE oligos as of August 2024 **(Oligo Design File)**
* 20240809_BARD1_SNVlib_ref_seqs_intron_annotated.xlsx - Reference sequences for each BARD1 SNV library with exonic and intronic sequence annotated **(BARD1 Lib Reference File)**
* 20241018_Orthogonal_BARD1_FunctionalAssays.xlsx - Curated list of variant previously assayed in orthogonal assays. Variants are listed along with the assay performed and performance of the variant **(BARD1 Orthogonal Assays File)**
* 20250128_BARD1_ClinVar_SNVsOnly_1Starplus.txt - Current BARD1 SNVs with at least one star **(BARD1 ClinVar file)**
* 20240905_BARD1_gnomADv4.1.0_SNVs.xlsx - BARD1 SNVs in gnomAD v4.1.0 **(BARD1 gnomAD file)**
* 20250123_BARD1_NormReadDepth.tsv - Normalized read depth across all BARD1 SGE targets for each replicate **(BARD1 Depth File)**
* 20250205_deletion_inputs.xlsx - BARD1 deletiion inputs with exons annoted for each region **(Depth Input File)**
* 20240830_BRCA1_SGE_AllScores.xlsx - BRCA1 SGE scores from Findlay et. al 2018 **(BRCA1 SGE File)**

Also contained are sub-directories used for specific analyses. These sub-directories are:
* ATG_Lib_data - Contains files associatd with the BARD1_X1A double start codon mutant experiment
  * Counts folder - contains TSV files that have counts for variants seen in the X1A double start experiment. (NC, D5/13 DNA, and RNA) **(ATG counts folder)**
  * X1A_annotations.xlsx - annotations used to annotate final score file **(X1A annotations file)**
  * 20250409_BARD1_X1A_ATG_scored.xlsx - Final score output file for variants in the double start experiment **(X1A Start Score file)**
* SNV_filtering_inputs - Contains input files needed to run the FilterImpossibleAAs script
    * 20250724_bard1.fasta - Genomic sequence for *BARD1* in fasta format **(BARD1 fasta file)**
    * 20250415_BARD1_filter_entry.xlsx - Input file for generating reference file **(BARD1 input file)**
    * 20240820_CURRENT_SGEoligos_2023.xlsx - Contains all BARD1 SGE oligos **(Oligo design file)**
    * 20240809_BARD1_SNVlib_ref_seqs_intron_annotated.xlsx - Reference file to be used in FilterImpossibleAAs **(lib reference file)**


