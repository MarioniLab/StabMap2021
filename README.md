# StabMap2021
#### Scripts used to analyse and generate figure panels for 'StabMap: Mosaic single cell data integration using non-overlapping features'.

The StabMap R package can be found [here](https://github.com/MarioniLab/StabMap).

## Scripts

The `scripts` folder contains an initialisation script `initialise.R` and functions for analysis and visualisation.

## Analysis

The `example` folder contains Rmarkdown documents used to analyse and generate figure panels:

- PBMC Multiome example and comparison:
    - `PBMC_Multiome_example.Rmd`
- Mouse Gastrulation scRNA-seq:
    - `MGD_chimera_generate.Rmd`
    - `MGD_StabMap_example.Rmd`
- PBMC Multiome non-overlapping example:
    - `PBMC_Multiome_nonoverlapping_example.Rmd`
- Mouse Gastrulation chimera and seqFISH analysis:
    - `MGD_seqFISH_generate.Rmd`
    - `MGD_chimera_StabMap_SeqFISH_example.Rmd`
    - `MGD_chimera_StabMap_SeqFISH_example_downstream.Rmd`
- Misc:
    - `StabMap_workflow_figures.R`
    - `UINMF_vignette_comparison.Rmd`
    - `MultiMAP_vignette_comparison.Rmd`

These documents assume the existence of folders outside of the working directory, `../Figures/raw/` for graph outputs and `../output/` for analysis output files.

## Contact

shila.ghazanfar \<at\> cruk.cam.ac.uk or marioni \<at\> ebi.ac.uk

