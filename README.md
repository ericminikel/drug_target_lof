This repository holds the data and code for the following manuscript:

[Minikel EV, Karczewski KJ, Martin HC, Cummings BB, Whiffin N, Alfoldi J, Trembath RC, van Heel DA, Daly MJ, Schreiber SL, MacArthur DG. Evaluating potential drug targets through human loss-of-function genetic variation. bioRxiv 530881. 2019 Jan 28.](https://www.biorxiv.org/content/10.1101/530881v1)

The R script [src/drug_target_lof_analysis.R](src/drug_target_lof_analysis.R) will reproduce the figures and statistics from the manuscript.

1. **System requirements:** We have run this code in R 3.5.1 and its dependencies are [`sqldf`](https://cran.r-project.org/web/packages/sqldf/index.html), [`plotrix`](https://cran.r-project.org/web/packages/plotrix/index.html), and the [exac_2015 repository](https://github.com/macarthur-lab/exac_2015).
2. **Installation guide:** Clone this repository, and edit lines 1 and 7 of the drug_target_lof_analysis.R script to point to the proper paths on your system.
3. **Running the code:** Simply run the drug_target_lof_analysis.R script to re-generate all of the figures. It will load all the necessary data from the [data](https://github.com/ericminikel/drug_target_lof/tree/master/data) directory in this repository itself. Statistics quoted in the manuscript text will appear in script output. The script runs in approximately 3 minutes on a MacBook Air.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

