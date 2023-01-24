# Veterinary Student Graduate Expectations
Code release supporting the analysis of survey data for a submitted manuscript. This repository includes data and code required to run a Shiny app that can be used to interactively conduct statistical explorations and tests from the survey results. Follow these steps to run the app:
1. Place `app.R` and `AnalysisDataSet.tar.gz` into a directory `/thisproject` 
2. Unzip `AnalysisDataSet.tar.gz` so that `AnalysisDataSet.csv` resides in `/thisproject`
3. Launch the Shiny app from R by running:
```
shiny::runApp('/thisproject')
```

Tested on R version 4.1.0. Requires packages 'nlme' and 'shiny'.

[![DOI](https://zenodo.org/badge/592933437.svg)](https://zenodo.org/badge/latestdoi/592933437)
