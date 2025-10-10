# BARD1_SGE_analysis

This repository contains all Python notebooks and scripts used to create figures, tables, and numerical figures used in the BARD1 SGE paper. All notebooks/scripts are organized into one of two folders. Notebooks and Scripts. Datasets used to generate figures are located in the Data folder. Unless otherwise stated, all notebooks and scripts used for figures can be made using the BARD1_final_table.xlsx found in the final_tables sub-folder in the Data folder.

Required supplmentary data tables and paths are already set for all notebooks. To regenerate figures, simply clone this repository. For figures generated through Python scripts, paths will need to be updated based on your local machine. 

## Notebooks
The Notebooks folder contains Python Notebooks that create individual panels used in the final figures and the supplmentary tables. The vizualization that will be created by each notebook is noted in the name of the notebook. Code in notebooks have been commented out and each contains a Markdown header describing the figure or table it produces

These notebooks are:
* SGE_AbnormalVCF - Generates .VCF file used as the input to Ensembl's VEP tool to get variant effect predictor annotations for AlphaMissense, REVEL, CADD, and SpliceAI **(Needs SGE SNV file)**
* SGE_ATGlib_analysis - Generates continuous log-fold change scores for variants in the BARD1_X1A double start codon mutant experiment as well as basic visualizations of the results **(Needs ATG counts folder and X1A annotations file)**
* SGE_ATGlib_HeatMap - Generates heat map visualizations of BARD1_X1A double start codon mutant experiment faceted by canonical start (Figure 5D)
* SGE_BARD1_DomainAnalysis - Generates strip plot showing distribution of missense scores across BARD1's domains (Figure 2C)
* SGE_BARD1vBRCA1 - Comparative analysis of missense sensitivity of BARD1 and BRCA1 domains **(Needs BRCA1 SGE File)**
* SGE_CaseControlAnalysis - Case-control analysis of BARD1 variants using the BRIDGES and CARRIERS breast cancer case-control cohorts. Output used in figure 3D and the BARD1_OddsRatios_table supplmentary table **(Needs BRIDGES and CARRIERS data)**
* SGE_ClinVar_analysis - Generates strip plot of ClinVar variants seen in SGE and does ROC analysis (Figures 3B and 3C)
* SGE_CorrelationMatrix_Heatmap - Generates Pearson's r boxplots across timepoints and regions and a Pearson r heatmap across all SGE sub-targets (Figure S1B and S1C).
* SGE_DeletionFigure_ReadDepth - Generates 3 figures: 1) Line plot of median read depth across all coding regions in BARD1 (figure S4), 2) 3bp deletion map across BARD1 (figure 2A) and 3) heatmap of scores across BARD1 (figure 2B)
* SGE_EditRate_BarPlot - Generates SGE target-facted bar plot showing editing rate generating useable reads across all targets and replicates (Figure S1A).
* SGE_Histogram_Stripplot - Creates histogram (top) and strip plot (bottom) of all BARD1 SGE scores (Figure 1B and 1C)
* SGE_MAFvsSGEScore - Creates 2D heatmap of MAF (retrieved from gnomAD and Regeneron Million Exomes) vs. SGE score (Figure 3A)
* **SGE_MakeFinalDataTable** - Builds BARD1_final_table supplementary table - the required input for all figure generating notebooks. Requires all files in /BARD1_SGE_analysis/Data/supp_table_inputs 
* SGE_OrthogonalAssays - Generates histogram and strip plot showing distribution of SGE scores for variants assayed in orthogonal functional assays (Figure S5B)
* SGE_PaperNumberGenerator - Generates numerical values used throughout the main text
* SGE_RNAanalysis - Builds scatter plot of RNA score vs. SGE score (Figure 5A) and RNA score vs. SGE score stem plots (Figure 5B, S8) 
* SGE_VEPs_vs_SGE - Creates scatter plots of predictor score vs. SGE score (Figure S5A) and case-control analysis of variants with moderate evidence towards benignity/pathogenicty from predictors (Figure S5C)
* SGE_ZnBindingAnalysis - Analysis of missense variants impacting Zn2+ coordinating residues in BARD1's RING domain. Generates strip plot (Figure S9) and comparative analysis with BRCA1 **(Needs BRCA1 SGE File)**
  
## Scripts
The Scripts folder contains script used to generate counts for the BARD1_X1A_ATG double mutant experiment and all scripts used to generate figures in PyMOL. Scripts are labled by the figure that will be generated and code is commented out.

These scripts are:
* BARD1_X1A_ATG_counter - Generates counts file for all variants in the BARD1_X1A_ATG double mutant experiment. 
* SGE_BARD1_colorPyMOL_MIS_only - Colors ribbon cartoon protein structures using missense SGE scores (Figure 6A)
* SGE_BRCA1_colorPyMOL_MIS - Colors ribbon cartoon protein structures for BRCA1 using missense SGE scores (Figure 6A) **(Needs BRCA1 SGE file)**
* SGE_DrawHelicalWheel - Generates helical-wheel diagrams for BARD1 and BRCA1 (Figure 6B)
* SGE_PyMOL_BARD1_ColorSurface - Generates colored space-filling models displaying minimum and mean missense surface score  (Figure 6D)
* SGE_PyMOL_BRCA1_ColorSurface - Generates colored space-filling models dispalying either minuimum or mean missense score at surface (Figure 6D) **(Needs BRCA1 SGE file)**
* SGE_PyMOL_Sphere - Generates ribbon cartoon protein structures with gray spheres at positions of interest where LoF missense variants exist (Figure 6C)


## Data
The data folder contains sub-directories containing all supplementary tables and data needed to generate figures. Data used to recompile provided supplementary data tables is also provided. 

## case_control_data
Contains data accessed from the BRIDGES and CARRIERS breast cancer case-control studies. Data for each study is organized into its respective folder:
* BRIDGES_data - Contains all data files from the BRIDGES study
  * 20250815_BRIDGES_missense_all.xlsx - All missense variants sequenced in the BRIDGES study
  * 20250815_BRIDGES_missense_population.xlsx - Missense variants only from population-based studies that were part of the BRIDGES study
  * 20250815_BRIDGES_PTVs_all.xlsx - All PTVs sequenced in the BRIDGES study
  * 20250815_BRIDGES_PTVs_pop.xlsx - PTVs only from population-based studies that were part of the BRIDGES study
 
* CARRIERS_data - Contains data file from the CARRIERS study
  * 20250303_CARRIERS_data.xlsx - Contains all variants sequenced in the CARRIERS case-control study.
 
## final_data_tables
Contains all final supplmentary tables and external data used for analyses and figures
* BRCA1_SGE_data.xlsx - Contains BRCA1 SGE data from Findlay et al. 2018 and Dace et al. 2025 in separate tabs
* BARD1_SGE_final_table.xlsx - Large supplmentary table containing all data from this study in multiple tabs. Used as input for all figure generation and analysis
* BARD1_OddsRatios_table.xlsx - Supplmentary table containing raw number of variants going into the case-control analysis done utilizing the BRIDGES and CARRIERS cohorts done in this study. 

### supp_table_inputs
Files used to create the BARD1_final_table supplementary table are found in the **supp_table_inputs** sub-directory. Files are dated based on when they were originally created/when data was accessed These files are: 
* 20240802_BARD1_Regeneron_MAF.xlsx - MAF data from Regeneron's Millon Exome study
* 20240905_BARD1_gnomADv4.1.0_SNVs.xlsx - MAF data for BARD1 SNVs accessed from gnomAD v4.1.0
* 20241016_Orthogonal_BARD1_FunctionalAssays.xlsx - Curated list of variant previously assayed in orthogonal assays. Variants are listed along with the assay performed and performance of the variant
* 20241217_BARD1_sgRNA_cutsites.xlsx - List of 1-based, hg38 genomic coordinates for Cas9 cut sites based on designed sgRNAs for each SGE target
* 20250415_BARD1_filter_entry.xlsx - Contains sheet of genomic coordinates for all SGE targets
* 20250825_BARD1counts.tsv - Contains raw counts of variants across all SGE targets and replicates
* 20250912_BARD1_ClinVarDels_1StarPlus.txt - Contains in-frame deletion variants accessed from Clinvar with at least 1-star review stuatus
* 20250912_BARD1_ClinVarSNVs_1StarPlus.txt - Contains SNVs accessed from Clinvar with at least 1-star review stuatus
* 20250922_BARD1RNAscores.tsv - Raw RNA scores for SNVs.
* 20250926_BARD1.editrates.sorted.tsv - Editing rates generating useable reads for each SGE target and replicate
* 20251002_BARD1delscores.tsv - Raw 3bp deletion scores
* 20251002_BARD1modelparam.tsv - Output parameters from GMM modeling. Gives best estimated thresholds for functional class classification
* 20251002_BARD1_snvs_VEP.xlsx - Ensembl VEP annotated SNVs files
* 20251002_BARD1snvscores.tsv - Raw SNV scores
* 20252002_BARD1snvscores.vcf - .VCF file input for Ensemble VEP annotation of SNVs (to variant effect predictor scores)
* 20251006_BARD1_MutPred2.xlsx - MutPred2 scores for coding missense SNVs
* 20251008_PhyloP.xlsx - PhyloP scores across BARD1

The supp_table_inputs sub-directory also dcontains two additional folders
* ATG_lib_data - Contains all data needed to recompute results from the X1A_ATG double mutant experiment
  * X1A_annotations.xlsx - Contains variant annotations for library X1A
  * 20250409_BARD1_X1A_ATG_scored.xlsx - Raw scores from the X1A_ATG double mutant experiment 
  * counts - Folder contains .tsv files with raw counts of variants for each replicate
* depth_data - Contains all data needed to recompute normalized depth values
  * deletion_inputs.xlsx - Contains coordinates for start and stop of each BARD1 coding exon
  * depth_files - Folder contains read depths at each base for each library and replicate




