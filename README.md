# StabMap2021
#### Scripts used to analyse and generate figure panels for 'StabMap: 
Stabilised mosaic single cell data integration using unshared features'.

The StabMap R package can be found [here](https://github.com/MarioniLab/StabMap).

## Scripts

The `scripts` folder contains an initialisation script `initialise.R` and 
functions for analysis and visualisation.

## Analysis

The `example` folder contains Rmarkdown documents used to analyse and 
generate figure panels:

- PBMC Multiome example and comparison:
    - `PBMC_Multiome_example.Rmd`
- Mouse Gastrulation scRNA-seq:
    - `MGD_chimera_generate.Rmd`
    - `MGD_StabMap_example.Rmd`
    - `MGD_Multihop_example.Rmd`
- PBMC Multiome non-overlapping example:
    - `PBMC_Multiome_nonoverlapping_example.Rmd`
- PBMC CyTOF, ECCITE-Seq and Multiome example:
    - `PBMC_CYTOF_generate.Rmd`
    - `PBMC_ECCITE_generate.Rmd`
    - `PBMC_Multihop_example.Rmd`
- Breast cancer IMC, CITE-Seq and Xenium example:
    - `BreastCancer_IMC_generate.Rmd`
    - `BreastCancer_CITE_generate.Rmd`
    - `BreastCancer_Xenium_generate.Rmd`
    - `BreastCancer_Multihop_example.Rmd`
- Mouse Gastrulation chimera and seqFISH analysis:
    - `MGD_seqFISH_generate.Rmd`
    - `MGD_chimera_StabMap_SeqFISH_example.Rmd`
    - `MGD_chimera_StabMap_SeqFISH_example_downstream.Rmd`
- Misc:
    - `StabMap_workflow_figures.R`
    - `UINMF_vignette_comparison.Rmd`
    - `MultiMAP_vignette_comparison.Rmd`

These documents assume the existence of folders outside of the `scripts` 
working directory, `../Figures/raw/` for graph outputs and `../output/` 
for analysis output files.

## Contact

shila.ghazanfar \<at\> cruk.cam.ac.uk, shila.ghazanfar \<at\> sydney.edu.au, or marioni \<at\> ebi.ac.uk

