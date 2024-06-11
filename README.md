# anemone-dynamics
Modelling body size changes in Nematostella and Aiptasia

Experimental measurements are in `Data/`, in a combination of Excel and CSV formats.

Linear and simple changepoint fits
====

The first set of code explores linear relationships and relatively simple dynamic models between the scientific variables. This is all performed in `main-analysis.R`, which calls `multiphasemodel.R` for helper functions modelling the multiphase dynamics. Most of the code is linear model fits and plots, but a bootstrapping section where several different experimental datasets are fit to several different changepoint model structures will take some time (perhaps half an hour on a single core at time of writing). The model fits are summarised and output to an automatically-generated Excel file; a summary figure of all multiphase model fits is also output by default. Other plots are generated on the fly for visualisation and validation.

Simulated annealing fits
====

The second set of code attempts to fit dynamic models to data using simulated annealing. This allows more complex changepoint models but also serves as validation for the simpler dynamics in the first set. This pipeline is wrapped by the Bash script `run.sh` which also serves to illustrate the workflow and intermediate files.

`export-data.R` and `aiptasia.R` pull original data into input formats for processing with the simulated annealing model fits. 

`optim-csv.c` attempts to fit changepoint models of various structures to input data, using simulated annealing to maximise likelihood.

`plot-results.R` outputs the results of these simulated annealing runs.
