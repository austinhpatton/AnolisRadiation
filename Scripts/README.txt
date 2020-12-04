This directory holds essential scripts. 

- Bootstrap-Contemp-Disparity.R runs the analysis to obtain estimates of contemporary disparity, and run the bootrapping analysis to assess significance (Figures 2C-F, S3).
- ConvergenceDivergence-Broad.R runs the analysis to assess the extent of convergence across the four major groups described in the manuscripts (SLA, GA, M1, M2: Figures 4F, S9-S12, Tables S2, S7, Supplementary Data 1 (in SI)).
- ConvergenceDivergence-Mainland.R runs a similar to analysis, but with a focus on mainland groups (Figure 6A).
- ConvergenceDivergenceFunctions.R - includes the functions sourced by the above two scripts to perform convergence analyses. 
- EcomorphologicalANCOVAs.R - performs phylogenetic ancovas (Figure 4A-E).
- HiSSE-Analysis.R - performs HiSSE analyses assuming empirical sampling fractions (Figures 2B, S18, Tables S5-S6).
- HiSSE-Analysis-Assume70perc.R - same as above, but assuming we have undersampled extant mainland diversity (Figure S19). 
- Mainland-ByRegion-Boxplots-ForM1-M2.R - runs the simple comparisons to test for differences among mainland regions in body size (Figure S15).
- Morphospace-FINAL.R - estimates hypervolumes (Figure 4A, 6C S7, S14).
- PhylogeneticPC-BAMM - summarized bamm runs using phylogenetic PCS (Figure 1A, S5). 
- PulledDR-Analysis.R - estimates the pulled diversification rates through time (Figure 3C-D).  
- TestAbundance-Conv-AlloSymp.R - Tests for the frequency of convergence among different biogeographic regions on the mainland, and performs a test of significance (Figure 6A). 
- Trait-Ecol-Boxplots.R - Tests for differences among major groups (GA, M1, M2) in morphology and ecology and visualizes as boxplots (Figure 5, S2, S8, S16, S17). 

The following are control files for the BAMM analyses:
Bamm-Diversification-ControlFile.txt - for estimating speciation rates.
BAMM-pPC_Control_File.txt - for estimating morphological rate of evolution using phylogenetic PCs. 
