[![DOI](https://zenodo.org/badge/346238990.svg)](https://zenodo.org/badge/latestdoi/346238990)

# yukon-chinook-diversity
Data and code to reproduce analyses in:
> Population diversity in Canadian-origin Chinook salmon: portfolio quantification and implications for conservation, harvest, equity, and resilience. In press. Connors B.M., Siegle M.R., Harding J., Rossi S., Staton B., Jones M., Bradford M., Browne R., Bechtol B., Doherty B. and S. Cox. Ecological Applications.

This project consists of 3 primary components: 
1. a state-space run-reconstruction model fitted to daily estimates of total border passage and population composition to estimate annual population-specific border passage and spawner abundance; 
2. age-structured state-space spawner-recruitment models to characterize diversity in population dynamics; and 
3. various analyses to quantify the consequences of population diversity for portfolio effects and equilibrium trade-offs between harvest and overfishing risk in the system.

The posterior samples from fitting the single- and multi-stock spawner-recruitment model are not included in this Git repo, as the file sizes are too large. Users may either fit the model themselves or request the file from the authors.

Project made possible by grant 1701 from the Arctic-Yukon-Kuskokwim Sustainable Salmon Initiative

## Folders and files
- `load.R`: Loads packages and scripts necessary for analysis. This file should be sourced prior to running other scripts.
- `01_inputs`: Data and model output to reproduce all figures in manuscripts.
- `02_run-reconstruction`: Code and data to fit the multi-stock run-reconstructions and generate some supplemental figures.
- `03_single-stock-SS-SRA`: Code to fit single stock state-space spawner-recruitment models. 
- `04_figures`: Code to generate manuscript main text figures. 
- `05_integrated-SS-SRA`: Code and data to fit a integrated population state-space spawner-recruit model and produce the content in Appendix S2 to the main-text.


