# DDE_BD
A R code and example data for paper titled "Bayesian inference of distributed time delay in transcriptional and translational regulation"

We propose a Bayesian inference framework based on replacing uninteresting or unobserved reactions with time delays. 
Although the resulting models are non-Markovian, recent results on stochastic systems with random delays allow us to 
rigorously obtain expressions for the likelihoods of model parameters. In turn, this allows us to extend MCMC methods 
to efficiently estimate reaction rates, and delay distribution parameters, from single-cell assays. We illustrate the advantages, 
and potential pitfalls, of the approach using a birth-death model with both synthetic and experimental data, and show that we can 
robustly infer model parameters using a relatively small number of measurements. We demonstrate how to do so even when only the 
relative molecule count within the cell is measured, as in the case of fluorescence microscopy.

This repository provides a R code and example data sets in order to estimate reaction constants of simple Birth-Death process with time delay. Gamma distribution is used for the delay time distribution and shape and rate parameter of gamma delay distribution are also estimated using this R code. Estimates of mean and variance of delay time are also presented.
