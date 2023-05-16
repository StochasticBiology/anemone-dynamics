# anemone-dynamics
Modelling body size changes in Nematostella and Aiptasia

Experimental measurements are in `Data/`, in a combination of Excel and CSV formats. `Precomputed/` stores pre-computed output from bootstrap resampling (see below). The analysis code is in two parts. 

Linear and simple changepoint fits
====

The first set of code explores linear relationships and relatively simple dynamic models between the scientific variables.

`Data_preparation_and_computations.ipynb` is a Jupyter notebook. The first part takes the original data and casts it into various input files for subsequent processing. Subsequent parts take the results from this processing (in the R files below) and analyse and visualise the results in various ways. 

`SimpleLinearRegression.R` explores the simple, one-phase relationships between variables, plotting visualisations and reporting inferred parameters and model selection statistics.

`MultiPhaseModels.R` fits changepoint-style models to the multi-phase experiments and performs bootstrapping for uncertainty quantification. This, because of the bootstrapping, takes a while. Its output also requires some manual work -- different experiments are analysed in the code, and the current version requires the filename for each experiment to be specified manually (see code comments). To address this, the Jupyter notebook by default reads precomputed output from the `Precomputed` directory.

Simulated annealing fits
====

The second set of code attempts to fit dynamic models to data using simulated annealing. This allows more complex changepoint models but also serves as validation for the simpler dynamics in the first set. This pipeline is wrapped by the Bash script `run.sh` which also serves to illustrate the workflow and intermediate files.

`export-data.R` and `aiptasia.R` pull original data into input formats for processing with the simulated annealing model fits. There is redundancy between these and the Jupyter notebook.

`optim-csv.c` attempts to fit changepoint models of various structures to input data, using simulated annealing to maximise likelihood.

`plot-results.R` outputs the results of these simulated annealing runs.
