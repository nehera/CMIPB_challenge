# CMIPB_challenge

INTRODUCTION
This repository includes our team's code used and submission for the 2nd CMI-PB prediction challenge: https://www.cmi-pb.org/blog/prediction-challenge-overview/

Team members included faculty and PhD students from the University of Minnesota's Biostatistics and Health Data Science Division:
* Aidan Neher
* Aparna Srinivaan
* Jamie Forschmiedt
* Aidan Dunleavy
* Steffen Ventz, PhD
* Eric Lock, PhD

We ultimately applied the [blockForest](https://bmcbioinformatics-biomedcentral-com.ezp2.lib.umn.edu/articles/10.1186/s12859-019-2942-y)<sup>1</sup> method to batch corrected training data and, using the normalized test data, predicted each challenge task's outcome. 

ABOUT THE FILES
* data_description.txt - this is a description of data features available for the prediction challenge. 
* missingness_EFL.R and missing_predict_EFL.R - these files reshape the processed data p.data provided by the challenge organizers and perform mean imputation for all features used in our model fitting. 
* random_forest.R - this file includes some experimentation with naive random forest and blockForest methods and includes 5-fold CV estimation of out-of-sample Spearman's correlation. 
* BlockRandomForest_EFL.r - this is what we used to train our blockForest models used in prediction of all test data tasks. 

REFERENCES
1. Hornung, R. & Wright, M. N. Block Forests: random forests for blocks of clinical and omics covariate data. BMC Bioinformatics 20, 1–17 (2019).
