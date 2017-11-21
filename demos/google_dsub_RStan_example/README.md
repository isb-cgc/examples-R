GooglePipelinesRStan

#################################################################################
# How to run a custom script with dsub
# David L Gibbs
# dgibbs@systemsbiology.org
# November 7, 2017

https://cloud.google.com/genomics/overview

#################################################################################

In this example, I'm going to be fitting Bayesian logistic regression models using Stan
(http://mc-stan.org/). Each job (model fitting) will process a single file,
but we could also have each job represent a parameter set, or model, all
processing the same data.

This collection of files can be used to demonstrate:

1. Generating a 'task matrix', each row describing a job in the google cloud.
   (cmd_generator.R, task_matrix.tsv)

2. Running a custom R script on user data.
   (stan_logistic_regression.R, data/*)

4. Using the Google dsub to automatically start up a VM, run a script, and shutdown.
Please see the 'how_to_dsub.txt' file for instructions.
