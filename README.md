here is the code for all steps of the analysis, divided into 

- OLINK_preprocessing_QC_PCA: checking the quality of the olink data ( repeated measurements,plates,missing rates,...), imputating NPX values under LOD, PCA plots, supplement figures (size/charge)

- ClinicalData_Prepartion: clinical data preparation: defintion of hypertension, extracting medication info at baseline, descriptive statistics

- CKDendpoint_definition_eGFRprogression: eGFR progression and definition of composite endpoint

- LASSO regression/Variance Explained Distributixon over Variable Category,for figure 1

- LinearModels_Contrasts.Rmd: figure 2, eGFR predictor for protein levels in linear models , and influence/differences in country and diagnose with contrasts 

- metadeconfoundR run: assocation of composite endpoint reached at different time points with olink proteins and confounding/assocaitons wth other clinical variables

- RandomForest.Rmd : Random Forest to predict reaching the ckd endpoint within a certain time frame (annual visits): data preparation, construct k-cv RFs and validation/visualization

- SurvivalAnalysis_Preparation.Rmd: defining survival time (time till endpoint is reached), calculating inflammation score and plot Kaplan-Meier curves for the i-scores as well as boxplots for egfr/alburie distribtution across i-scores

- CoxModesl.Rmd : cox regression model and visualization 

- protein_network_iscore_subgroups.Rmd: constructing protein correlation network, analysing clinical parameters in the inflammation score subgroups
