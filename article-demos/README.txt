==========================
This .zip archive contains the replication files for the JSS submission:

"hdMTD: An R Package for High-Dimensional Mixture Transition Distribution Models"
Authors: Maiara Gripp, Giulio Iacobelli, Guilherme Ost, Daniel Y. Takahashi
==========================

Contents:
- code.R                          : Main replication script.
- code.html                       : Precompiled HTML output (using knitr::spin()) for quick inspection.
- simulated_data.rds              : Simulated dataset used in examples, preprocessed for efficiency.
- results_sequential_selection.rds: Preprocessed results for a time-consuming analysis.
- hdMTD_outputs.rds		  : Preprocessed results for a time-consuming estimations.		

==========================
How to run the code
==========================

1. Open R and set the working directory to this folder.
2. Install required packages if not already installed:
   dplyr, ggplot2, lubridate, purrr, tidyr, hdMTD, future and future.apply.

3. Run the script:
   knitr::spin("code.R")

This will generate a .md and .html file.

====================================
Important notes on reproducibility:
====================================

Some lines in code.R are preceded by a logical variable "recompute <- FALSE". 
These correspond to computations that may take several minutes or more to complete.
For convenience, their results have been precomputed and saved in the .rds files
listed above under "Contents".

- For a quick inspection, you may simply leave all "recompute <- FALSE" lines
unchanged; in this case, the precomputed results will be loaded.
- To recompute a result simply set the corresponding line to "recompute <- TRUE".
- However, to fully reproduce all results from scratch, you may set the default
"recompute_all <- FALSE" at the beginning of the code to "recompute_all <- TRUE".
- All instructions are provided within code.R itself, along with approximate time durations
(on an i7-1255U 10-core processor) for these time-consuming chunks.
- To overwrite the RDS files set `save_precomputed <- TRUE` in the beginning of code.R.

====================================
About the included code.html file:
====================================

The file code.html in ReplicationMaterials.zip was generated with all computations fully executed — that is,
 "recompute_all <- TRUE". Therefore, no preprocessed data was used.


==========================
Contact
==========================

For questions, please contact:
Maiara Gripp  
<maiara@dme.ufrj.br>