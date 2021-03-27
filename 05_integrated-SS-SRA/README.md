This subdirectory contains the necessary files to fit an alternative approach to the mixed-stock spawner-recruitment analysis we performed for the analysis presented within the manuscript. Details of both approaches we used can be found in Supplement B of the manuscript.

The posterior samples from fitting this model are not included in this Git repo, as the file size is too large. Users may either fit the model themselves ([JAGS 4.3.0](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/) used) or request the file from the authors.

The directory structure is as follows:

1. `fit-integrated-model`: this subdirectory contains the data, functions, and other code to fit the integrated  population spawner-recruitment analysis (originally presented in Staton et al. [2020](https://cdnsciencepub.com/doi/abs/10.1139/cjfas-2019-0281)) to Canadian-origin Chinook salmon.
2. `sra-supplement`: this subdirectory contains the code and supporting files to create Supplement B to the main text. Supplement B contains model descriptions for both the spawner-recruitment analysis in the main text and for the integrated approach, as well as detailed model output and comparisons between these two approaches.

