==========================
This .zip archive contains the replication files for the JSS submission:

"hdMTD: An R Package for High-Dimensional Mixture Transition Distribution Models"
Authors: Maiara Gripp, Guilherme Ost, Giulio Iacobelli, Daniel Y. Takahashi
==========================

Contents:
- code.R                          : Main replication script.
- code.html                       : Precompiled HTML output (using knitr::spin()) for quick inspection.
- simulated_data.rds              : Simulated dataset used in examples, preprocessed for efficiency.
- results_sequential_selection.rds : Preprocessed results for a time-consuming analysis.

==========================
How to run the code
==========================

1. Open R and set the working directory to this folder.
2. Install required packages if not already installed:
   dplyr, ggplot2, lubridate, purrr, tidyr, hdMTD

3. Run the script:
   knitr::spin("code.R")

This will generate a .md and .html file.

====================================
Important notes on reproducibility:
====================================

Some lines in code.R are commented out and marked with a "**WARNING**" label. 
These correspond to computations that may take several minutes or more to complete.

- These lines can be ignored for a quick inspection.
- However, to fully reproduce all results from scratch, 
  you may manually uncomment and run these sections.
- All instructions are provided within code.R itself, along with approximate time durations.

Additionally, some analyses use preprocessed results stored in the .rds files for efficiency.
These cases are explicitly indicated in code.R, and instructions for regenerating the results
from raw simulations are also included, along with approximate runtimes on i7-1255U, 10 cores.

====================================
About the included code.html file:
====================================

The file code.html in ReplicationMaterials.zip was generated with all computations fully executed — that is, no sections were commented out, and no preprocessed data was used.

==========================
Contact
==========================

For questions, please contact:
Maiara Gripp  
<maiara@dme.ufrj.br>