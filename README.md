# DDE_BD

Although the resulting models are non-Markovian, recent results on stochastic systems with random delays allow us to 
rigorously obtain expressions for the likelihoods of model parameters. In turn, this allows us to extend MCMC methods 
to efficiently estimate reaction rates, and delay distribution parameters, from single-cell assays. We illustrate the advantages, 
and potential pitfalls, of the approach using a birth-death model with both synthetic and experimental data, and show that we can 
robustly infer model parameters using a relatively small number of measurements. We demonstrate how to do so even when only the 
relative molecule count within the cell is measured, as in the case of fluorescence microscopy.

This repository provides the toolkit for a statistical inference mathod in order to estimate reaction constants of simple Birth-Death process with time delay for the paper on ?????????. Gamma distribution is used for the delay time distribution and inferecne for shape and rate parameter of gamma delay distribution are also included. Estimates of mean and variance of delay time are presented as well as reaction constants and parameter of gamma delay distribution.

File description 
Birth.Death.Delay.R: R code for a statistical inference for Birth-Death process with delay
YFP1.rda: R example data for experiment data set 1 (40 trajectories)
YFP2.rda: R example data for experiment data set 2 (29 trajectories)
YFP1.rda: R example data truncated at T=22 of YFP1
YFP1.rda: R example data truncated at T=22 of YFP2
