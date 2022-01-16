# Digital-Data-Analysis

From both folder Altoida_VAMBN_paper and ADNI_VAMBN_paper
Prepare the data for VAMBN model:
This is the main analysis file. 
  1. run the R files clean_format->impute_aux->data_types (scripts with fixed settings)
  2. Grid search to generate best hyperparameters
  3. run HI-VAE jupyter notebook (up until VP decoding)
  4. run full R script bnet.R bayesian network, bn prediction evaluation of DMs and MMSE (only for Altoida) and VirtualPatient validation
  5. decode generated VPs in the HI-VAE notebook
  6. run R scritpt finalRandomForestPred.R for random forest prediction evaluation of DMs and MMSe (only for Altoida)

Evaluation of Virtual Patients
  1. Run R script plot_all_marginals_paper.R
  2. Run R script correlations.R and correlations_paper.R

For Classifiers
  1. Run R scripts classifierSGL.R and classifierSGLind.R
  2. Comparisons of performance of classifiers aucAll.R and aucCLassifier.ipynb

For scatter plot generation between different features
  1. Run scatterPlots.R

Extra for ADNI
  1. dataPreparationADNI.Rmd to process ADNI data before clean_format.R
  1. To generate DMs in ADNI -> counterfactual_cognitiveDomain.R and counterfactual_cognitiveDomain.R before bnet.R




