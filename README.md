# typhoidAMR

How to Run

1)	To run the model simulations across different R0 estimates, proportion symptomatic, vaccination coverage, and initial resistance, run outputScriptSim.m. This will generate SampleResPCs.mat, which contains the regression model resDiff_mdl used to estimate the change in proportion resistant across 10 years of vaccine roll-out, and DiffVec, a vector of these differences.
a.	The default is resistance emerging 10 years before vaccine roll-out begins. To generate the 5 and 20 year scenarios, edit line 42 in sir_amr_age.mat
2)	Next, run OutputTablesDR_ResDiff.m. This will import several data files: 
a.	file containing estimates for FQNS and MDR in each of the 73 countries included in this study
b.	file containing vaccine coverage for each country
c.	R0, proportion symptomatic estimates for each country from Bilcke et al. 2019
d.	SampleResPCs.mat
The file will generate DRtable.mat, which contains mean, min and max estimates for FQNS and MDR for each country (region when not available); and resDiffstable.csv, which contains 2000 sample estimates for each country of change in proportion resistant.
3)	resDiffsTable.csv is then used along with the beta distribution fits for FQNS and MDR (FQNDists.csv, MDRDists.csv) in a series of R scripts to generate estimates of baseline and Camp15 cases, deaths and DALYs.
4)	To generate these tables, user runs TCEA_INPUT. With no edits, this will generate baseline and Camp15 estimates for total typhoid. User can generate estimates for FQNS, MDR, or Total AMR by commenting out appropriate lines (38-43) in RAWDATA_amrams.r, and commenting in the lines that source the _resonly.r files in TCEA_INPUT and commenting out the original lines (lines 97 for 89, 167 for 164, 184 for 179)
5)	This will generate files called summary_baseline10yrs_undisc.csv and summary_RC15nodisc_table.csv. These should be renamed with the appropriate file name based on the scenario, e.g. if MDR only was run, the files are renamed Baseline10Years_73_MDR.csv and Camp15Impact_73_MDR.csv. There are files with rounded estimates for ease of reading, and unrounded estimates for analysis (with suffix NR). 
6)	These files are transferred back to matlab, and imported in in TablesFigs_BetaDist.m.
7)	TablesFigs_BetaDist.m uses these estimates to generate  Figures 1 and 2 in the main text, and Figures S5-9.

