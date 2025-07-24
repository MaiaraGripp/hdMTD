\# Replication Files for JSS Submission



\*\*"hdMTD: An R Package for High-Dimensional Mixture Transition Distribution Models"\*\*  

Authors: Maiara Gripp, Guilherme Ost, Giulio Iacobelli, Daniel Y. Takahashi



---



\## Contents



\- `code.R`: Main replication script  

\- `code.html`: Precompiled HTML output (using `knitr::spin()`) for quick inspection  

\- `simulated\_data.rds`: Simulated dataset used in examples, preprocessed for efficiency  

\- `results\_sequential\_selection.rds`: Preprocessed results for time-consuming analysis  



---



\## How to Run the Code



1\. Open R and set the working directory to this folder.  

2\. Install required packages if not already installed: 

   dplyr, ggplot2, lubridate, purrr, tidyr, hdMTD



3\. Run the script:



knitr::spin("code.R")



This will generate a .md and .html file.



\# IMPORTANT NOTES:



Some lines in code.R are commented out and marked with a "\*\*WARNING\*\*" label.

These correspond to computations that may take several minutes or more to complete.



\- These lines can be ignored for a quick inspection.

\- However, to fully reproduce all results from scratch,

  you may manually uncomment and run these sections.

\- All instructions are provided within code.R itself along with approximate time durations.

\- Some analyses use preprocessed results stored in the .rds files for efficiency.

\- It is explictly informed in code.R when this happens and the instructions to regenerate these files from raw simulations are also included in code.R along with approximate time duration.





