# PPRforEstimatingAutoregressiveOperator
Creating an R package for estimating autoregressive operator using Principal Predictive Components (PPC) according to Kargin-Onatski algorithm.

# Running the code
Simply go to the R script `Test.R`, change the working directory and run it as usual.
It contains a simply test to check that the package has been installed, and then a comparison between the package and an R function for the PPC.
Data are taken from `CanadianWeather`. Although the error is yet quite big, performance are improved yet for MSE yet for computational time.
This big MSE could be due to the small amount of time instants?

# What's missing
The code has been developed only for the algorithm given a parameter (CV version not yet implemented).
