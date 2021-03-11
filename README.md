# yukon-chinook-diversity
> Population diversity in Canadian-origin Chinook salmon: portfolio quantification and implications for conservation, harvest, equity, and resilience. In preparation. Connors B.M., Siegle M.R., Harding J., Rossi S., Staton B., Jones M., Bradford M., Browne R., Bechtol B., Doherty B. and S. Cox.

Data and code to:
1. estimate annual population-specific border passage and abundance via a state-space run-reconstruction model fitted to daily estimates of total border passage and population composition; 
2. characterize diversity in population dynamics by fitting age-structured, state-space spawner-recruitment models; and 
3. quantify the consequences of population diversity for portfolio effects and equilibrium trade-offs between harvest and overfishing risk

Project made possible by grant 1701 from the Arctic-Yukon-Kuskokwim Sustainable Salmon Initiative

## Folders and files
- `01_inputs`: Data neccessaryu to reprouce all figures in manuscripts
- `figures.R`: Generate figures in manuscript
- 
-  `load.R`: Loads packages and scripts necessary for analysis. This file should be sourced prior to running other scripts.

- `functions.R`: All functions written for the analysis should be placed in this file.
  
- `close_loop_sims.R`: Run closed loop forwad simulations.


  
- `simulation_summary.Rmd`: R Markdown doc that summarizes simulations, and generates figures for manuscript (in `images` folder)

- `appendix_A.R`: R Markdown doc that details the stationary Ricker to time-varying Beverton-Holt formulation we used, and simulations to justify its parameterization in our closed-loop simulations. 
