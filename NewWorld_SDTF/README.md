# NewWorld_SDTF
R scripts to conduct analyses in "The phylogenetic diversity and structure of the seasonally dry forests in the Neotropics" by Arango et al. Submitted to Journal of Biogeography

File description: 

**MasterScript.R**: this script allows to follow the steps used in the analyses and calculations performed in this study. It is separated in sections by #### and it has comments on most important actions using #. The packages used in each analysis are called at the start of each section.


**allsites.csv**: This file contains all the sites dowloaded from the DRYFLOR data base, with a special column indicating whether it was used or not in our analyses. Only the sites used in our study have the "code" and "bioregion" columns filled.

**biogeo_codes.csv**: Cointains all the sites used in our analysis with an unique numeric code ("codigo") used as ID. 

**Biogeo_sp.csv**: All the woody species present in the DRYFLOR inventories used in our study. The "A_N column" represents the accepted names (by 2016) of the species according to APG IV.

**ALLMBMOD.txt**: A modified version of the Smith & Brown (2016) megaphylogeny with Magallón (2015) and tree of life backbones. The modifications were made to some species that were misspelled. References cited in the main text.
