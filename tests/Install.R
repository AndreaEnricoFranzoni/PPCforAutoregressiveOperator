#to uploda the packages
#change here the directory
#setwd("/Users/andreafranzoni/Documents/Politecnico/Magistrale/PACS/PPCforAutoregressiveOperator")
#then 
#Rcpp::compileAttributes(".") 
# and then push the changes to modify the exports

library(devtools)
devtools::install_github("AndreaEnricoFranzoni/PPCforAutoregressiveOperator", force = TRUE)
